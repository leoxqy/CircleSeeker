# Module Renames

Python module names use underscores. This table maps legacy/internal names to the new, descriptive names (with their Chinese descriptions). Where hyphenated names are shown, the corresponding Python import uses underscores.

- carousel -> tandem-to-ring (`circleseeker.modules.tandem_to_ring`): 串联重复序列环化分析
- gatekeeper -> um-classify (`circleseeker.modules.um_classify`): U/M分类筛选
- trapeze -> cecc-build (`circleseeker.modules.cecc_build`): 复杂eccDNA重建
- menagerie -> umc-process (`circleseeker.modules.umc_process`): U/M/C整合处理
- harmonizer -> ecc-dedup (`circleseeker.modules.ecc_dedup`): 基于CD-HIT的去重复
- sieve -> read-filter (`circleseeker.modules.read_filter`): 读段过滤
- contortionist -> split-refine (`circleseeker.modules.split_refine`): 分裂读段区域精炼
- juggler -> circle-validate (`circleseeker.modules.circle_validate`): 环状DNA验证
- playbill -> report-generator (`circleseeker.modules.report_generator`): HTML报告生成
- propmaster -> file-organizer (`circleseeker.modules.file_organizer`): 文件整理组织

Notes:
- Old imports continue to work. The new names are aliases re-exporting the same objects.
- If CLI tools with hyphenated names are desired, we can add console scripts to `pyproject.toml` that wrap the corresponding functionality.

