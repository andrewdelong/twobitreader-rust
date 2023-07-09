// Standard library
use std::iter::once;
use std::mem::transmute;
use std::panic;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;

// Stack size for spawned threads.
// Small stack suffices for usage in this crate, but not for general use.
const WORKER_STACK_SIZE: usize = 32 * 1024;

// Calls a closure on each element of an iterator, in parallel.
//
// If any invocation of f() panics, the initial thread will also panic.
//
// This function is designed for use within the create, and has not been
// designed or tested for more general use cases. The goal is simple
// std-based data parallelism that foregoes dependency on rayon.
//
pub fn parallel_for<I, F, T>(iter: I, f: F)
where
    I: ExactSizeIterator<Item = T> + Sync + Send, // ExactSize makes bounding num_threads simple
    F: Fn(T) + Sync + Send,
{
    // Define a dummy error type so that try_parallel_for can be used as the implementation.
    #[derive(Debug)]
    struct UnusedError;

    // Reuse try_parallel_for but ignore the result since there's no way to return an error.
    try_parallel_for(iter, |item| -> Result<(), UnusedError> {
        f(item);
        Ok(())
    })
    .unwrap();
}

// Calls a closure on each element of an iterator, in parallel.
//
// If any invocation of f returns an error, the result of this function
// will contain the first error detected.
//
// If any invocation of f panics, the calling thread will also panic.
//
// This function is designed for use within this create, and has not been
// designed or tested for more general use cases. The goal is simple
// std-based data parallelism, without dependency on rayon.
//
// TODO: When #![feature(try_trait_v2)] becomes stable, design this
// around the Try trait rather than Result.
//
pub fn try_parallel_for<I, F, T, E>(iter: I, f: F) -> Result<(), E>
where
    I: ExactSizeIterator<Item = T> + Sync + Send, // ExactSize makes bounding num_threads simple
    F: Fn(T) -> Result<(), E> + Sync + Send,
    E: Sync + Send + 'static,
{
    // Worker loop. Consumes the next item from mutex's iter, then applies f(item) to it.
    let consume_iter = |mutex: Arc<Mutex<I>>| -> Result<(), E> {
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
    let num_threads = num_cpus::get().min(iter.len()).max(1);
    let mutex = Arc::new(Mutex::new(iter));
    let join_handles = (0..num_threads - 1)
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

    // Call consume_iter once from current thread (i.e., reuse it as a worker), then join
    // and collect all Result values. Exit with 'panic code' (101) if any thread panicked.
    let join_results = once(Ok(consume_iter(mutex)))
        .chain(join_handles.into_iter().map(|join_handle| join_handle.join()))
        .collect::<Vec<_>>();

    // Check that all thread joins succeeded, and extract the worker results (of consume_iter).
    let worker_results = join_results
        .into_iter()
        .map(|join_result| match join_result {
            Ok(worker_result) => worker_result,
            Err(panic_error) => panic::resume_unwind(panic_error),
        })
        .collect::<Vec<_>>();

    // Check that all workers succeeded. If not, return the first error encountered.
    match worker_results.into_iter().find_map(|r| r.err()) {
        Some(worker_error) => Err(worker_error),
        None => Ok(()),
    }
}
