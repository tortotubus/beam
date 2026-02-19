typedef int sig_atomic_t; // TODO: add to defaults.h

static volatile sig_atomic_t g_shutdown_requested = 0;
static volatile sig_atomic_t g_shutdown_signal = 0;

static void shutdown_handler (int signo) {
  g_shutdown_requested = 1;
  g_shutdown_signal = (sig_atomic_t) signo;
}

int install_shutdown_handlers (void) {
  struct sigaction sa;
  memset (&sa, 0, sizeof (sa));
  sa.sa_handler = shutdown_handler;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;

  // const int sigs = checkpoint_configuration.checkpoint_on_signals;
  for (size_t i = 0; i < sizeof (checkpoint_configuration.checkpoint_on_signals) / sizeof (checkpoint_configuration.checkpoint_on_signals[0]); ++i) {
    if (sigaction (checkpoint_configuration.checkpoint_on_signals[i], &sa, NULL) != 0) {
      return -1; // errno is set
    }
  }
  return 0;
}

int shutdown_requested (void) {
  return (int) g_shutdown_requested;
}

int shutdown_signal (void) {
  return (int) g_shutdown_signal;
}
