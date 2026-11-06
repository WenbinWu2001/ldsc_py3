"""Help formatting shared by CLI entry points and workflow parsers.

Keep option names and paths intact when descriptions wrap so users can read
and copy dependencies without reconstructing tokens split across lines.
Shared path descriptions distinguish chromosome suites from ordinary globs;
workflow parsers use only the descriptions supported by their input readers.
The options every materializing command shares, ``--output-dir`` and
``--log-level``, are added by one helper each, so their type, choices and
default are stated once; each command supplies its own help text.
"""

import argparse
import textwrap

from .workflow_logging import LOG_LEVEL_CHOICES, LOG_LEVEL_HELP


CHROMOSOME_PATH_HELP = (
    "Quote patterns: '*' matches filename text; '@' substitutes chromosome numbers 1-22 by default. "
    "Use '*' to select available files for a chromosome subset where allowed."
)
SCALAR_PATH_HELP = (
    "Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded."
)



def output_directory_argument(value: str) -> str:
    """Parse a destination directory option, rejecting a blank value that would mean the working directory."""
    if not value.strip():
        raise argparse.ArgumentTypeError(
            "must name a directory; an empty value would write into the working directory "
            "(pass '.' to use the working directory on purpose)"
        )
    return value


def add_output_dir_argument(group, *, help: str) -> None:
    """Add the required ``--output-dir DIR`` option, which rejects a blank directory, with ``help``."""
    group.add_argument('--output-dir', required=True, metavar='DIR', type=output_directory_argument, help=help)


def add_log_level_argument(group, *, help: str = LOG_LEVEL_HELP) -> None:
    """Add the ``--log-level`` option (default ``INFO``) with the shared level choices."""
    group.add_argument('--log-level', default='INFO', choices=LOG_LEVEL_CHOICES, help=help)


class CLIHelpFormatter(argparse.HelpFormatter):
    """Wrap help at spaces without splitting option names or file paths."""

    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(" ", text).strip()
        return textwrap.wrap(text, width, break_long_words=False, break_on_hyphens=False)

    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(" ", text).strip()
        return textwrap.fill(text, width, initial_indent=indent, subsequent_indent=indent,
                             break_long_words=False, break_on_hyphens=False)
