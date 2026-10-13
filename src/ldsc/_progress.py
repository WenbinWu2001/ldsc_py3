"""Low-volume, elapsed-time progress records for workflow phases.

Progress stays on the existing LDSC logger routes. Call once per bounded chunk
or object; periodic records are time-throttled and never inspect input data.
"""

from time import perf_counter
from collections import Counter
from contextvars import ContextVar
from functools import wraps
import logging


class PhaseProgress:
    """Log phase boundaries and periodic completed counts with constant memory."""

    def __init__(self, logger, phase, label, total=None, *, interval=30.0, level=logging.INFO):
        self.logger, self.phase, self.label = logger, phase, label
        self.total, self.interval = total, interval
        self.completed = 0
        self.object = ''
        self.level = level
        self.outcomes = Counter()

    def _log(self, event):
        counts = f'; completed={self.completed}' if self.completed or self.total is not None else ''
        if self.total is not None:
            counts += f'/{self.total}'
        outcomes = ''.join(f'; {name}={count}' for name, count in sorted(self.outcomes.items()))
        current = f'; object={self.object}' if self.object and self.logger.isEnabledFor(logging.DEBUG) else ''
        self.logger.log(self.level, 'Phase %s: %s; %s%s%s%s; elapsed=%.2fs.',
                         event, self.phase, self.label, counts, outcomes, current, perf_counter() - self.started,
                         extra={'phase': self.phase, 'phase_event': event, 'phase_label': self.label,
                                'completed': self.completed, 'total': self.total, 'current_object': self.object,
                                'outcomes': dict(self.outcomes)})

    def __enter__(self):
        self.started = self.last = perf_counter()
        self._log('start')
        return self

    def advance(self, count=1, *, object='', force=False, outcome=None):
        self.completed += count
        if outcome is not None:
            self.outcomes[outcome] += count
        self.object = object
        now = perf_counter()
        if force or now - self.last >= self.interval:
            self._log('progress')
            self.last = now

    def __exit__(self, kind, value, traceback):
        self._log('failed' if kind else 'complete')
        return False


_ACTIVE_PHASE = ContextVar('ldsc_progress_phase', default=None)


def report_phase(logger, phase, label, *, level=logging.INFO):
    """Add phase boundaries to a synchronous workflow without changing its return."""
    def decorate(function):
        @wraps(function)
        def reported(*args, **kwargs):
            with PhaseProgress(logger, phase, label, level=level) as progress:
                token = _ACTIVE_PHASE.set(progress)
                try:
                    return function(*args, **kwargs)
                finally:
                    _ACTIVE_PHASE.reset(token)
        return reported
    return decorate


def advance(count=1, *, object='', force=False, total=None):
    """Report bounded-chunk progress when a caller has installed a phase."""
    progress = _ACTIVE_PHASE.get()
    if progress is not None:
        if total is not None:
            progress.total = total
        progress.advance(count, object=object, force=force)
