"""Standard exit codes for CircleSeeker CLI.

Following shell conventions:
- 0: Success
- 1: General/business logic error
- 2: Command line usage error
- 130: Terminated by SIGINT (128 + 2)
- 143: Terminated by SIGTERM (128 + 15)
"""

EXIT_SUCCESS = 0
EXIT_ERROR = 1  # General/business logic error
EXIT_USAGE = 2  # Command line usage error
EXIT_SIGINT = 130  # 128 + SIGINT(2)
EXIT_SIGTERM = 143  # 128 + SIGTERM(15)
