// #include <unistd.h>
#include <errno.h>
#include <stdio.h>
#include <sys/stat.h> // mkdir
#include <sys/types.h>

bool file_exists (const char* fname) {
  FILE* fp = fopen (fname, "rb");
  if (!fp) {
    return false;
  } else {
    fclose (fp);
    return true;
  }
}

int copy_file (const char* src, const char* dst) {
  FILE* in = fopen (src, "rb");
  if (!in) {
    perror ("fopen src");
    return -1;
  }

  FILE* out = fopen (dst, "wb");
  if (!out) {
    perror ("fopen dst");
    fclose (in);
    return -1;
  }

  char buf[1 << 16]; // 64 KB
  size_t n;

  while ((n = fread (buf, 1, sizeof (buf), in)) > 0) {
    if (fwrite (buf, 1, n, out) != n) {
      perror ("fwrite");
      fclose (in);
      fclose (out);
      return -1;
    }
  }

  if (ferror (in)) {
    perror ("fread");
    fclose (in);
    fclose (out);
    return -1;
  }

  fclose (in);
  fclose (out);
  return 0;
}

int create_path (const char* path) {
  if (mkdir (path, 0777) == 0)
    return 0;
  if (errno == EEXIST)
    return 0;
  return -1;
}

static char executable_name[4096] = "";

const char* get_executable_name (void) {
  if (executable_name[0] == '\0') {
    ssize_t n = readlink (
      "/proc/self/exe", executable_name, sizeof (executable_name) - 1);
    if (n < 0) {
      // fallback
      return "unknown";
    }
    executable_name[n] = '\0';
  }

  const char* slash = strrchr (executable_name, '/');
  return slash ? slash + 1 : executable_name;
}