
```sh
isoforms <- detect_isoform("C:/Users/bc/data/isomiR/meta.csv", "vvi")
isomirs <- generate_isomir(isoforms)
expr <- generate_expr(isomirs)
tsi <- get_tsi(expr)

mirnas <- load_mirna("vvi", 13)
sam <- to_sam(isoforms, mirnas)
```
