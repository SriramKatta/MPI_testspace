import pandas as pd
import matplotlib.pyplot as plt
import sys

# Read CSV
df = pd.read_csv(sys.argv[1])

# Clean column names
df.columns = df.columns.str.strip().str.lstrip('#').str.replace(' ', '_')

# Now group by problem_size
for size, group in df.groupby("problem_size"):
    group = group.sort_values("procs")
    plt.plot(group["procs"], group["flop_rate"], marker="o", label=f"Size {size}")

plt.xlabel("Processors count")
plt.ylabel("FLOP Rate(Mfolps/sec)")
plt.title("Strong Scaling FLOP Rate vs Number of Processors")
plt.legend(title="Problem Size")
plt.grid(True)

# Start x-axis at 0
# plt.xlim(left=0)

# Only place ticks at actual processor numbers
proc_ticks = sorted(df["procs"].unique())
plt.xticks(proc_ticks)

plt.tight_layout()
plt.savefig(f"{sys.argv[1]}.png", dpi=300)
