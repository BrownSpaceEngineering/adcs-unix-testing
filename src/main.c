#include "include/test.h"
#include <stdio.h>

int main(void) {
    int failures = test_run_all();
    // Non-zero exit status on failure so scripts/CI can tell, but every test still runs first
    return failures == 0 ? 0 : 1;
}
