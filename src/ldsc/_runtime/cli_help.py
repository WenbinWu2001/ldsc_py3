"""Help formatting shared by CLI entry points and workflow parsers.

Keep option names and paths intact when descriptions wrap so users can read
and copy dependencies without reconstructing tokens split across lines.
Shared path descriptions distinguish chromosome suites from ordinary globs;
workflow parsers use only the descriptions supported by their input readers.
"""

import argparse
import textwrap


CHROMOSOME_PATH_HELP = (
    "Quote patterns: '*' matches filename text; '@' substitutes chromosome numbers 1-22 by default. "
    "Use '*' to select available files for a chromosome subset where allowed."
)
SCALAR_PATH_HELP = (
    "Accepts an exact path or a quoted '*' pattern matching exactly one file; '@' is not expanded."
)


class CLIHelpFormatter(argparse.HelpFormatter):
    """Wrap help at spaces without splitting option names or file paths."""

    def _split_lines(self, text, width):
        text = self._whitespace_matcher.sub(" ", text).strip()
        return textwrap.wrap(text, width, break_long_words=False, break_on_hyphens=False)

    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(" ", text).strip()
        return textwrap.fill(text, width, initial_indent=indent, subsequent_indent=indent,
                             break_long_words=False, break_on_hyphens=False)
