echo "#ifndef COMMIT_H\n#define COMMIT_H\n\n    const std::string Commit = \""$(cat .git/refs/heads/master)\" ";\n\n#endif // COMMIT_H" > 2L_SDVRP/include/Commit.h
