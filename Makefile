# Compiler flags
# Include git branch and commit hash in the build
CFLAGS += -D'GIT_BRANCH_NAME="$(GIT_BRANCH_NAME)"' -D'GIT_COMMIT_HASH="$(GIT_COMMIT_HASH)"'

# Robust error checking
CFLAGS += -Wextra -Werror -Werror=maybe-uninitialized
CFLAGS += -Wshadow -Wnull-dereference -Wduplicated-cond -Wlogical-op -Werror=return-type -Wfloat-equal
CFLAGS += -Wdangling-else -Wtautological-compare
CFLAGS += -fwrapv # Enable fwrapv (wrap on overflow of signed integers) just to be safe
CFLAGS += -fsigned-char # Ensure that char is signed as your average c programmer might expect -- it's actually default unsigned on arm!

# Linking to Ccontrol 
# CFLAGS += -L$(PATH_TO_CCONTROL) -lCControl

# Disable warnings for unused parameters due to ASF functions having unused parameters
CFLAGS += -Wno-unused-parameter #Because some ASF functions have unused parameters, supress this warning

INCLUDE_DIRS = -I./ -I./linalg/LinearAlgebra/ -I./linalg/Lapack/Include/
CFLAGS += $(INCLUDE_DIRS)

.phony = all clean 

EXECS = adcs_test 
DEPS = test.c 
HEADS = test.h 

all: $(EXECS)

adcs_test: $(HEADS) $(DEPS)
	gcc $(CFLAGS) $^ -o $@

clean: 
	rm -f $(EXECS)