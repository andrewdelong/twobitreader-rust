// Standard library
use std::mem::transmute;
use std::panic;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;

// Stack size for worker threads.
// Modest stack suffices for usage in this crate, but not for general use.
const WORKER_STACK_SIZE: usize = 32 * 1024;

// Calls f() on each item of an iterator in parallel.
//
// If any invocation of f() returns an error, the result of this function
// will contain the first error detected.
//
// If any invocation of f() panics, the panic is propagated to the calling thread.
//
// This function is designed for the simple use case within this crate, and is not
// optimized or tested for other use cases. The only reason for this code is a simple
// std-based parallel loop, without dependency on rayon.
//
// TODO: Consider using Try trait when #![feature(try_trait_v2)] is stable.
//
pub(crate) fn try_parallel_for<I, F, E>(iter: I, f: F) -> Result<(), E>
where
    I: ExactSizeIterator + Sync + Send,
    F: Fn(I::Item) -> Result<(), E> + Sync + Send,
    E: Sync + Send + 'static,
{
    // Worker loop. Consumes the next item from mutex's iter, then applies f(item) to it.
    let consume_iter = |mutex: Arc<Mutex<I>>| {
        loop {
            // Consume next item from mutex-protected iter, then release the lock with drop()
            let mut iter = mutex.lock().unwrap();
            let item = iter.next();
            drop(iter);

            // Apply f(item) or terminate if no more items.
            match item {
                Some(item) => f(item)?,
                None => break Ok(()),
            }
        }
    };

    // Spawn additional threads to consume items in parallel.
    let num_threads = thread::available_parallelism().unwrap().get().min(iter.len()).max(1);
    let mutex = Arc::new(Mutex::new(iter));
    let join_handles = (0..num_threads)
        .map(|_| {
            // Clone reference to iter mutex
            let mutex = Arc::clone(&mutex);

            // Override static lifetime requirement on spawn's closure.
            // Doing so is ONLY SAFE if ALL THREADS are joined BEFORE RETURNING from this function.
            // Until #![feature(thread_spawn_unchecked)] moves out of nightly and into stable,
            // this transmute hack from crossbeam will have to do.
            let closure: Box<dyn FnOnce() -> Result<(), E> + Send> = Box::new(move || consume_iter(mutex));
            let closure: Box<dyn FnOnce() -> Result<(), E> + Send + 'static> = unsafe { transmute(closure) };

            // Spawn using a modest stack size. Plenty for the use cases within this crate.
            thread::Builder::new().stack_size(WORKER_STACK_SIZE).spawn(closure).unwrap()
        })
        .collect::<Vec<_>>();

    // Collect all thread join results.
    let join_results = join_handles.into_iter().map(|h| h.join()).collect::<Vec<_>>();

    // Check that all thread joins succeeded. If not, propagate the first panic.
    if join_results.iter().any(|r| r.is_err()) {
        panic::resume_unwind(join_results.into_iter().find_map(|r| r.err()).unwrap());
    }

    // Check that all worker results (after unwrapping join results) are Ok. If not, return the first error.
    if let Some(worker_error) = join_results.into_iter().find_map(|r| r.unwrap().err()) {
        return Err(worker_error);
    }

    Ok(())
}
