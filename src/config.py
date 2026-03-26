from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent
data_dir = Path("/work3/s252608/DL_project/data")

raw_data = data_dir / "raw"
processed_data = data_dir / "processed"

HVG_count = 5000