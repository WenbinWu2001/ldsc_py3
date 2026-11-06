"""One row of an input-issue repair table (``diagnostics/input_issues.tsv``).

Every input gate, including PLINK discovery in ``paths``, builds its rows
with :func:`input_issue`. The repair is a required argument, so each row
names the repair for its own input instead of inheriting another input's.
"""

INPUT_ISSUE_COLUMNS = ['input_role', 'source', 'chrom', 'reason', 'details', 'repair']


def input_issue(role, source, chrom, reason, details, *, repair):
    """Return one input-issue row in the ``INPUT_ISSUE_COLUMNS`` schema.

    ``source`` (a path or label) and ``details`` (text or an exception) are
    stored as text; ``chrom`` is the chromosome label, or ``''`` when the
    issue concerns no single chromosome.
    """
    return dict(input_role=role, source=str(source), chrom=chrom, reason=reason, details=str(details),
                repair=repair)
