library(ggplot2)
library(dplyr)
suppressPackageStartupMessages(library("argparse"))  
parser <- ArgumentParser()
parser$add_argument("--dir", help="directory of the job")
parser$add_argument("--jobid", help="Job Id")
args <- parser$parse_args()


# Cargar los datos desde el archivo
t <- read.table(paste(args$dir,'AuxFiles/DF_Colores', sep=""), sep = ',', header = TRUE)

# Crear el dataframe
df <- data.frame(
  PDBID = t$PDBID,
  pos = t$pos,
  AA = t$AA,
  FrstState= t$FstState,
  FstI=t$FstI
)

# Definir los colores
colores <- c("NEU" = "gray", "MAX" = "red", "MIN" = "green", "-" = "white")

# Filtrar las filas con valores NA en la columna FstI
df_filtered <- df[!is.na(df$FstI), ]

filas <- list()

# Obtener los PDBID únicos
pdbids <- unique(df$PDBID)

# Iterar sobre cada PDBID
for (pdbid in pdbids) {
  # Filtrar los datos por PDBID
  datos_pdbid <- df[df$PDBID == pdbid, ]
  
  # Crear una fila con los valores de FstI
  fila <- datos_pdbid$FstI
  
  # Agregar la fila a la lista
  filas[[pdbid]] <- fila
}

# Crear un nuevo data frame con las filas
df_resultado <- as.data.frame(do.call(rbind, filas))
rownames(df_resultado) <- pdbids

df_resultado[df_resultado == "N/A"] <- NA
df_resultado[is.na(df_resultado)] <- 0

# Realizar el clustering utilizando K-means
num_clusters <- as.integer(length(pdbids)/3)  # Número de clusters deseado
set.seed(123)  # Fijar semilla para reproducibilidad
clusters <- kmeans(df_resultado, centers = num_clusters)

# Obtener los resultados del clustering
df_resultado$cluster <- clusters$cluster

# Ordenar la lista por el clustering
df_ordenado <- df_resultado[order(df_resultado$cluster), ]

# Agregar columna de clustering al dataframe
df_filtered$clustering <- df_resultado$cluster
cluster_labels <- as.character(clusters$cluster)
df_filtered$PDBID <- factor(df_filtered$PDBID, levels = rev(unique(df_filtered$PDBID)))
lfilas = length(filas)
laa=length(t$AA)
if ( lfilas < 20){
  laa=laa*lfilas
  lfilas = 20
}
# Ordenar el gráfico por clustering
p <- ggplot(df_filtered, aes(x = pos, y = reorder(PDBID, PDBID), fill = FrstState)) +
  geom_tile(colour = NA, linewidth = 1) +
  geom_text(aes(label = AA), size = lfilas*0.03, color = "black") +
  scale_fill_manual(values = colores) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = lfilas*0.1),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 5),  # Ajusta el tamaño del título de la leyenda aquí
        legend.text = element_text(size = 5),
        legend.margin = margin(t = 1, unit = "lines"),
        legend.key.size = unit(0.2, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"))
warnings()
# Ajustar el tamaño de la leyenda en la parte superior
p <- p + theme(legend.box.margin = margin(t = -lfilas*0.5, unit = "pt"))  # Ajusta el tamaño de la leyenda aquí
h=lfilas*0.08
w=laa*0.0025
ggsave(paste(args$dir,'OutPutFiles/MSFA_',args$jobid,'.png', sep=""), p, width = w, height = h, dpi = 300, bg = "white",limitsize = FALSE)
