from common import TABLE_DIR
from Table_4 import main


target = TABLE_DIR / "Table_S3.csv"
if not target.exists():
    main()
else:
    print(f"Table S3 already exists: {target}")
