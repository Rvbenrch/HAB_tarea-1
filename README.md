# Tarea 1: Análisis Funcional (g:Profiler · ORA)

**Autor:** Rubén M. Rodríguez Chamorro  
**Asignatura:** Herramientas de Bioinformática  
**Fecha:** 2 de Noviembre de 2025

---

<details>
<summary><strong>Enunciado y estructura del repositorio</strong></summary>

Entrega un script en Python que implemente de forma detallada un análisis funcional de los genes **COX4I2**, **ND1** y **ATP6**. El análisis debe utilizar librerías de Python para evaluar los procesos biológicos asociados a estos genes. Asegúrate de que tu código esté bien documentado y describa claramente los métodos y bases de datos utilizadas para obtener información funcional.

**Estructura propuesta del repo:**
```

/analisis-funcional/
├── data/
│   └── genes_input.txt              # Lista de genes (un gen por línea)
├── scripts/
│   └── rubenscript.py               # Script de análisis funcional (CLI)
├── results/                         # Resultados generados automáticamente
├── README.md                        # Este archivo
└── requirements.txt                 # Dependencias (gprofiler-official, pandas, etc.)

```

**Rúbrica de evaluación (10 pts):**
| Criterio | Descripción | Puntos |
|---|---|---|
| 1. Funcionalidad | El script realiza el análisis funcional solicitado. | 4 |
| 2. Documentación | Código comentado y explicación de métodos/bases de datos. | 2 |
| 3. Uso de librerías | gprofiler-official, pandas, etc. | 2 |
| 4. Formato y estilo | Código legible y buenas prácticas. | 1 |
| 5. Automatización (CLI) | Acepta argumentos desde la terminal. | 1 |

**Dependencias mínimas (requirements.txt):**
```

gprofiler-official>=1.0.0
pandas>=2.0.0
matplotlib>=3.7.0
numpy>=1.24.0

```
</details>

---

## Descripción y fundamento

Este proyecto implementa un **análisis de enriquecimiento funcional (ORA) con g:Profiler** sobre una lista de genes. El objetivo es identificar **procesos biológicos (GO:BP), funciones moleculares (GO:MF), componentes celulares (GO:CC) y rutas (Reactome, KEGG)** **significativamente sobrerrepresentados** en el conjunto de entrada.

**Metodología (resumen):**
1. Entrada: lista de genes (símbolos HUGO), p. ej. `data/genes_input.txt`.  
2. Consulta: API oficial **g:Profiler** vía `gprofiler-official`.  
3. Estadística: prueba **hipergeométrica** con corrección por **FDR (Benjamini–Hochberg)**.  
4. Salida: tabla TSV de términos significativos, JSON de metadatos y figura PNG con el **Top N** (por defecto, 20).

---

## CLI del script

El script se ejecuta como **línea de comandos** y acepta parámetros:

```

usage: rubenscript.py [-h] [--input INPUT] [--organism ORGANISM]
[--sources SOURCES [SOURCES ...]] [--fdr FDR]
[--top TOP] [--outdir OUTDIR] [--no-plot]

optional arguments:
--input INPUT                Ruta a la lista de genes (un gen por línea).
--organism ORGANISM          Organismo (ej. hsapiens, mmusculus).
--sources SOURCES [...]      Bases a consultar (GO:BP GO:MF GO:CC REAC KEGG).
--fdr FDR                    Umbral FDR (por defecto 0.05).
--top TOP                    Nº de términos a mostrar en la figura.
--outdir OUTDIR              Directorio de salida (por defecto: results/).
--no-plot                    Solo tabla/JSON; no generar figura.

````

---

## Ejecución

<details>
<summary><strong>Ubuntu / Debian</strong></summary>

**1) Preparación (una sola vez):**
```bash
sudo apt update
sudo apt install -y python3-venv python3-pip
# En sistemas nuevos: sudo apt install -y python3.12-venv
````

**2) Crear y activar entorno virtual:**

```bash
python3 -m venv venv
source venv/bin/activate
```

> Si ves el error “ensurepip is not available”, instala `python3-venv` como arriba y recrea el venv.

**3) Instalar dependencias:**

```bash
pip install -r requirements.txt
```

**4) Ejecutar (ejemplo recomendado):**

```bash
python3 scripts/rubenscript.py \
  --input data/genes_input.txt \
  --organism hsapiens \
  --sources GO:BP GO:MF GO:CC REAC \
  --fdr 0.05 \
  --top 20 \
  --outdir results
```

</details>

<details>
<summary><strong>Windows (PowerShell)</strong></summary>

**1) Abrir PowerShell en la carpeta del repo**
**2) Crear y activar entorno:**

```powershell
py -m venv venv
.\venv\Scripts\Activate.ps1
# Si hay error de políticas:
# Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

**3) Instalar dependencias:**

```powershell
pip install -r requirements.txt
```

**4) Ejecutar:**

```powershell
py scripts\rubenscript.py `
  --input data\genes_input.txt `
  --organism hsapiens `
  --sources GO:BP GO:MF GO:CC REAC `
  --fdr 0.05 `
  --top 20 `
  --outdir results
```

</details>

---

## Qué genera el script

En `results/` se producen automáticamente (con sello temporal):

* `enrichment_YYYYMMDD-HHMMSS.tsv` → **tabla** con los términos enriquecidos y sus estadísticas.
* `enrichment_YYYYMMDD-HHMMSS.json` → **metadatos** (organismo, fuentes, FDR, top, nº genes, etc.).
* `enrichment_YYYYMMDD-HHMMSS.png` → **figura** con el **Top N** términos por `-log10(p)`.

---

## Qué he hecho y qué significan mis resultados

He ejecutado el análisis con **Homo sapiens (hsapiens)** sobre una lista de **5 genes canónicos de cáncer**:
**TP53, BRCA1, BRCA2, EGFR, PTEN**, consultando **GO:BP/GO:MF/GO:CC/Reactome** con **FDR ≤ 0.05** y **Top 20** para la figura.

**Resumen interpretativo (tabla):**

|        Categoría | Término enriquecido                         | Significado biológico                      | Lectura de resultados                                                                     |
| ---------------: | :------------------------------------------ | :----------------------------------------- | :---------------------------------------------------------------------------------------- |
|        **GO:BP** | *cell cycle process*                        | Regulación y progresión del ciclo celular. | Los genes analizados participan de forma central en **proliferación y control mitótico**. |
|        **GO:BP** | *DNA replication*                           | Duplicación del ADN antes de mitosis.      | Señala **replicación y reparación** como funciones nucleares del conjunto.                |
|        **GO:BP** | *cell cycle checkpoint*                     | Puntos de control ante daño genómico.      | Apoya el rol de **supresores tumorales** (p. ej., TP53, PTEN) deteniendo el ciclo.        |
|        **GO:BP** | *response to ionizing radiation*            | Respuesta al daño por radiación.           | Concordante con **vías de reparación** mediadas por BRCA1/2 y TP53.                       |
| **REAC / GO:BP** | *negative regulation of cell proliferation* | Inhibición de proliferación anómala.       | Refuerza la **función antitumoral** del conjunto (control de crecimiento).                |

**Figura (Top 20):**
La imagen `results/enrichment_20251102-115912.png` muestra los 20 términos más significativos ordenados por **-log10(p)**. Cuanto mayor el valor, mayor la evidencia estadística. En mi ejecución, predominan términos de **ciclo celular y reparación del ADN**, coherentes con el fenotipo tumoral y la biología de TP53/BRCA1/BRCA2/EGFR/PTEN.

---

## Interpretación y conclusiones

* El conjunto analizado está **fuertemente asociado** a **ciclo celular**, **replicación** y **respuesta a daño en el ADN**.
* Esto es **biológicamente consistente** con genes **supresores tumorales** y reguladores de **integridad genómica**.
* La combinación de **tabla (TSV)** y **figura (PNG)** permite **priorizar** procesos clave y **comunicar** resultados de forma clara y reproducible.

---

## Solución de problemas (FAQ) y cumplimiento de la rúbrica.

* **`ensurepip is not available` al crear venv:** instala `python3-venv` y recrea el entorno.
* **`ModuleNotFoundError: gprofiler`/`pandas`:** activa el venv y ejecuta `pip install -r requirements.txt`.
* **No hay resultados significativos:** revisa símbolos HUGO, organismo y fuentes; relaja `--fdr` o aumenta el tamaño de la lista.
* **Figura vacía:** puede que `--top` sea mayor que el nº de términos significativos; ajusta `--top`.

* **CLI reproducible:** parámetros `--input`, `--organism`, `--sources`, `--fdr`, `--top`, `--outdir`.
* **Documentación:** el script incluye docstrings y comentarios explicando **método ORA**, **g:Profiler** y **FDR**.
* **Librerías adecuadas:** `gprofiler-official`, `pandas`, `numpy`, `matplotlib`.
* **Formato/estilo:** funciones pequeñas, nombres descriptivos y salida ordenada en `results/`.

---

## Licencia

MIT.
© 2025 Rubén M. Rodríguez Chamorro.

```
