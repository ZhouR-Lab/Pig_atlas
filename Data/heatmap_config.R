# heatmap_config.R
# 存储热图绘制所需的固定顺序和颜色映射参数

# 1. 定义列的固定排序
specified_order <- c(
  "PN28#Kidney#Erythroid","PN180#Kidney#Erythroid","E38#Kidney#Erythroid","E80#Kidney#Erythroid","PN0#Kidney#Erythroid",
  "E38#Lung#Erythroid","E80#Lung#Erythroid","PN0#Lung#Erythroid","PN28#Lung#Erythroid","PN180#Lung#Erythroid",
  "E38#Liver#Erythroid","E80#Liver#Erythroid","PN0#Liver#Erythroid","PN28#Liver#Erythroid","PN180#Liver#Erythroid",
  "E38#Kidney#B cells","E80#Kidney#B cells","PN0#Kidney#B cells","PN28#Kidney#B cells","PN180#Kidney#B cells",
  "E38#Lung#B cells","E80#Lung#B cells","PN0#Lung#B cells","PN28#Lung#B cells","PN180#Lung#B cells",
  "E38#Liver#B cells","E80#Liver#B cells","PN0#Liver#B cells","PN28#Liver#B cells","PN180#Liver#B cells",
  "E80#Lung#Plasma","PN0#Lung#Plasma","PN28#Lung#Plasma","PN180#Lung#Plasma",
  "E38#Liver#NK/NKT","PN28#Kidney#NK/NKT","PN180#Kidney#NK/NKT","E80#Kidney#NK/NKT","PN0#Kidney#NK/NKT",
  "E80#Lung#NK/NKT","PN0#Lung#NK/NKT","PN28#Lung#NK/NKT","PN180#Lung#NK/NKT",
  "E80#Liver#NK/NKT","PN0#Liver#NK/NKT","PN28#Liver#NK/NKT","PN180#Liver#NK/NKT",
  "PN0#Kidney#T cells","PN28#Kidney#T cells","PN180#Kidney#T cells",
  "PN0#Lung#T cells","PN28#Lung#T cells","PN180#Lung#T cells",
  "PN0#Liver#T cells","PN28#Liver#T cells","PN180#Liver#T cells",
  "E80#Kidney#T cells","E80#Lung#T cells","E80#Liver#T cells","E38#Kidney#T cells","E38#Lung#T cells","E38#Liver#T cells",
  "E38#Lung#Mast cells","PN0#Lung#Mast cells","PN28#Lung#Mast cells","PN180#Lung#Mast cells",
  "E38#Liver#Megakaryocytes","E80#Liver#Megakaryocytes","PN0#Liver#Megakaryocytes","PN28#Liver#Megakaryocytes","PN180#Liver#Megakaryocytes",
  "E38#Kidney#DC","E80#Kidney#DC","PN0#Kidney#DC","PN28#Kidney#DC","PN180#Kidney#DC",
  "E38#Lung#DC","E80#Lung#DC","PN0#Lung#DC","PN28#Lung#DC","PN180#Lung#DC",
  "E38#Liver#DC","E80#Liver#DC","PN0#Liver#DC","PN28#Liver#DC","PN180#Liver#DC",
  "E38#Lung#Neutrophils","E80#Lung#Neutrophils","PN0#Lung#Neutrophils","PN28#Lung#Neutrophils","PN180#Lung#Neutrophils",
  "E38#Liver#Neutrophils","E80#Liver#Neutrophils","PN0#Liver#Neutrophils","PN28#Liver#Neutrophils","PN180#Liver#Neutrophils",
  "E80#Kidney#Monocytes","PN0#Kidney#Monocytes","PN28#Kidney#Monocytes","PN180#Kidney#Monocytes",
  "E80#Lung#Monocytes","PN0#Lung#Monocytes","E38#Liver#Monocytes","E80#Liver#Monocytes","PN0#Liver#Monocytes",
  "PN28#Liver#Monocytes","PN180#Liver#Monocytes","E38#Kidney#Monocytes","PN28#Lung#Monocytes","PN180#Lung#Monocytes",
  "E38#Liver#Kupffer","E80#Liver#Kupffer","PN0#Liver#Kupffer","PN28#Liver#Kupffer","PN180#Liver#Kupffer",
  "E38#Heart#Macrophages","E80#Heart#Macrophages","PN0#Heart#Macrophages","PN28#Heart#Macrophages","PN180#Heart#Macrophages",
  "E38#Kidney#Macrophages","E80#Kidney#Macrophages","PN0#Kidney#Macrophages","PN28#Kidney#Macrophages","PN180#Kidney#Macrophages",
  "E38#Lung#Macrophages","E80#Lung#Macrophages","PN0#Lung#Macrophages","PN28#Lung#Macrophages","PN180#Lung#Macrophages",
  "E38#Muscle#Macrophages","E80#Muscle#Macrophages","PN0#Muscle#Macrophages","PN28#Muscle#Macrophages","PN180#Muscle#Macrophages"
)

# 2. 定义注释的颜色映射字典
col_tissue <- c("Muscle"="#3498CB", "Lung"="#895b8a", "Heart"="#EB8C86", "Kidney"="#EED314", "Liver"="#339848")
col_stage  <- c("E38"="#430155", "E80"="#31688E", "PN0"="#20918D", "PN28"="#35B879", "PN180"="#FEE724")
col_cell   <- c("Erythroid"="#bf242a", "B cells"="#e0ebaf", "Plasma"="#c7dc68", "T cells"="#badcad",
                "NK/NKT"="#88cb7f", "Mast cells"="#bbbcde", "DC"="#E1BEE7", "Kupffer"="#f2a7da",
                "Neutrophils"="#9e2977", "Monocytes"="#BA68C8", "Macrophages"="#9C27B0", "Megakaryocytes"="#b44c97")