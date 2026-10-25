# Simulación y análisis de señales con la Transformada de Fourier

Proyecto en Python para generar y analizar un pulso rectangular, una función escalón y una senoide de 50 Hz. El programa calcula la FFT, grafica la magnitud y la fase del espectro, y comprueba propiedades de linealidad, desplazamiento temporal y escalamiento en frecuencia.

## Requisitos

- Python 3.9 o posterior.
- NumPy.
- Matplotlib.

## Instalar dependencias en Windows

Abre PowerShell o la terminal de Visual Studio Code en la carpeta del proyecto y ejecuta:

```powershell
py -m pip install -r requirements.txt
```

Si el comando `py` no está disponible, prueba:

```powershell
python -m pip install -r requirements.txt
```

## Ejecutar

Desde la carpeta donde está `transformada_fourier.py`, ejecuta:

```powershell
py transformada_fourier.py
```

Las gráficas se guardan automáticamente en una carpeta `resultados` junto al script. Para elegir otra carpeta de salida:

```powershell
py transformada_fourier.py --salida mis_graficas
```

## Archivos generados

- `resultados/01_senales_tiempo.png`: pulso, escalón y senoide en el dominio del tiempo.
- `resultados/02_espectros_magnitud_fase.png`: magnitud y fase de los espectros.
- `resultados/03_propiedades_fourier.png`: pruebas de linealidad, desplazamiento y escalamiento.

## Parámetros de simulación

Se usan 2048 muestras a una frecuencia de muestreo de 1024 Hz, en una ventana de 2 segundos. La resolución frecuencial es 0.5 Hz y la frecuencia de Nyquist es 512 Hz. La senoide tiene una frecuencia de 50 Hz; el pulso rectangular tiene amplitud 1 y ancho aproximado de 0.40 s.

El programa calcula una FFT de un solo lado para señales reales y escala el resultado por el intervalo entre muestras. La fase de componentes con magnitud casi nula se oculta porque no es estable. Como la FFT procesa un número finito de muestras, los límites de la ventana afectan el espectro; esto es especialmente visible para la función escalón.

## Publicar en GitHub

Incluye en el repositorio `transformada_fourier.py`, `README.md`, `requirements.txt` y la carpeta `resultados` con las gráficas. No es necesario subir el entorno virtual `.venv` ni la carpeta `__pycache__`.
