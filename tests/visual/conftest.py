"""Force a non-interactive backend before any visual test imports pyplot."""

import matplotlib

matplotlib.use("Agg")
