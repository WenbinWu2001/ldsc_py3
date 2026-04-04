# For pldsc, break annotation files input format into two sets: query annot and baseline annot. Each result is written separately.
# - query annot allows multiple annotations, each tested separately against baseline annot
# - baseline annot are read once and recycled to save comp time
# Add safe handling when pldsc fails for one query annot and continues to the next query annot, instead of crashing the whole program.