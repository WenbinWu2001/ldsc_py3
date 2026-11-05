"""Parse ``ldsc`` command lines with the parsers the CLI uses, without running them."""

from ldsc import cli
from ldsc._runtime.workflow_kwargs import parse_command_line


def parse_command(argv):
    """Parse ``argv`` (subcommand first) as ``ldsc`` does and record the subcommand as ``args.command``."""
    command, *rest = argv
    args = parse_command_line(cli._command_parser(command), rest)
    args.command = command
    return args
