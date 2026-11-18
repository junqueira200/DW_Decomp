#!/bin/sh
# Writes the current git commit hash to 2L_SDVRP/include/Commit.h.
# printf (not echo) so "\n" is a newline on both macOS and Linux.
commit=$(git rev-parse HEAD)
printf '#ifndef COMMIT_H\n#define COMMIT_H\n\n#include <string>\nconst std::string Commit = "%s";\n\n#endif // COMMIT_H\n' "$commit" > 2L_SDVRP/include/Commit.h
