/** @file elff.cpp
 * @brief Dispatch XML cases to the installed serial or MPI Basilisk drivers.
 */
#include <elff/io/C/xmlCase.h>

#include <cerrno>
#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

int fail(const std::string &message) {
  std::fprintf(stderr, "elff: %s\n", message.c_str());
  return 1;
}

// Inspect launcher metadata before exec 
bool launched_parallel() {
  for (const char *name : {"OMPI_COMM_WORLD_RANK", "PMI_RANK", "PMIX_RANK"}) {
    const char *value = std::getenv(name);
    int rank = -1;
    if (!value)
      continue;
    const char *end = value + std::strlen(value);
    const auto result = std::from_chars(value, end, rank);
    if (result.ec == std::errc{} && result.ptr == end && rank >= 0)
      return true;
  }
  const char *step = std::getenv("SLURM_STEP_ID");
  const char *value = std::getenv("SLURM_PROCID");
  if (!step || !*step || !value || !std::strcmp(step, "batch") ||
      !std::strcmp(step, "extern"))
    return false;
  int rank = -1;
  const char *end = value + std::strlen(value);
  const auto result = std::from_chars(value, end, rank);
  return result.ec == std::errc{} && result.ptr == end && rank >= 0;
}

void usage() {
  std::puts("Usage: elff case.xml [driver arguments...]\n"
            "The XML Basilisk.caseType selects the case; serial/MPI is "
            "detected automatically.\n"
            "Launch parallel runs with mpiexec or srun.\n"
            "  --help  Show this help.");
}

} // namespace

int main(int argc, char **argv) {
  try {
    std::string selected;
    int index = 1;
    if (argc > 1 &&
        (!std::strcmp(argv[1], "--help") || !std::strcmp(argv[1], "-h"))) {
      usage();
      return 0;
    }
    if (argc > 1 && !std::strcmp(argv[1], "--"))
      ++index;
    else if (argc > 1 && argv[1][0] == '-')
      return fail(std::string("Unknown option: ") + argv[1]);
    if (index == argc) {
      usage();
      return fail("An XML case file is required.");
    }
    const bool mpi = launched_parallel();
    elff_xml_document_t raw = nullptr;
    elff_xml_error error{};
    auto status = elff_xml_document_load_file(argv[index], &raw, &error);
    std::unique_ptr<elff_xml_document_handle,
                    decltype(&elff_xml_document_destroy)>
        document(raw, elff_xml_document_destroy);
    if (status != ELFF_XML_OK)
      return fail(error.message);
    elff_xml_node_t root = nullptr, basilisk = nullptr;
    status = elff_xml_document_root(raw, "Case", &root, &error);
    if (status != ELFF_XML_OK)
      return fail(error.message);
    status = elff_xml_node_unique_child(root, "Basilisk", &basilisk, &error);
    if (status == ELFF_XML_OK) {
      char name[64] = {};
      status = elff_xml_node_attribute_string(basilisk, "caseType", name,
                                              sizeof(name), nullptr, &error);
      if (status != ELFF_XML_OK)
        return fail(error.message);
      selected = name;
    } else if (status != ELFF_XML_NOT_FOUND)
      return fail(error.message);
    if (selected.empty())
      return fail("The dispatcher requires <Basilisk caseType=\"NAME\"> in the XML.");
    if (selected != "flapping_fibre" && selected != "settling_fibre")
      return fail("Unsupported XML case type: " + selected);

    const auto directory = std::filesystem::read_symlink("/proc/self/exe").parent_path();
    const auto executable = directory / (selected + (mpi ? "_mpi" : "_serial"));

    if (access(executable.c_str(), X_OK) != 0)
      return fail("Cannot execute " + executable.string() + ": " +
                  std::strerror(errno));
    std::vector<std::string> arguments;
    arguments.push_back(executable.string());
    for (; index < argc; ++index)
      arguments.emplace_back(argv[index]);
    std::vector<char *> command;
    for (auto &argument : arguments)
      command.push_back(argument.data());
    command.push_back(nullptr);
    document.reset();
    execvp(command[0], command.data());
    return fail("Cannot launch " + arguments.front() + ": " +
                std::strerror(errno));
  } catch (const std::exception &error) {
    return fail(error.what());
  }
}
