"""Parse ``ldsc`` command lines with the parsers the CLI uses, without running them."""

from ldsc import cli


def parse_command(argv):
    """Parse ``argv`` (subcommand first) and record the subcommand as ``args.command``."""
    command, *rest = argv
    args = cli._command_parser(command).parse_args(rest)
    args.command = command
    return args
