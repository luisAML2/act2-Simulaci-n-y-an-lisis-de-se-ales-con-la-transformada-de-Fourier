from pathlib import Path
import argparse

import matplotlib.pyplot as plt
import numpy as np



FS = 1024.0
N = 2048
DT = 1.0 / FS
T = (np.arange(N) - N // 2) * DT


def crear_senales(t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pulso = ((t >= -0.20) & (t <= 0.20)).astype(float)
    escalon = (t >= 0.0).astype(float)
    senoide = np.sin(2.0 * np.pi * 50.0 * t)
    return pulso, escalon, senoide


def transformar(x: np.ndarray, t: np.ndarray = T) -> tuple[np.ndarray, np.ndarray]:
    frecuencia = np.fft.rfftfreq(len(x), d=DT)
    espectro = np.fft.rfft(x) * DT
    espectro *= np.exp(-2j * np.pi * frecuencia * t[0])
    return frecuencia, espectro


def graficar_senales(pulso: np.ndarray, escalon: np.ndarray,
                     senoide: np.ndarray, carpeta: Path) -> None:
    fig, ejes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    configuracion = [
        (pulso, "Pulso rectangular", "Amplitud"),
        (escalon, "Función escalón", "Amplitud"),
        (senoide, "Función senoidal (50 Hz)", "Amplitud"),
    ]
    for ax, (x, titulo, ylabel) in zip(ejes, configuracion):
        ax.plot(T, x, lw=1.5)
        ax.set_title(titulo, loc="left", fontsize=11)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
    ejes[-1].set_xlabel("Tiempo (s)")
    ejes[-1].set_xlim(T[0], T[-1])
    fig.suptitle("Señales en el dominio del tiempo", fontsize=15)
    fig.tight_layout()
    fig.savefig(carpeta / "01_senales_tiempo.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def graficar_espectros(senales: list[tuple[str, np.ndarray]], carpeta: Path) -> None:
    """Guarda magnitud y fase de las tres señales."""
    fig, ejes = plt.subplots(3, 2, figsize=(12, 10), sharex="col")
    for fila, (nombre, x) in enumerate(senales):
        frecuencia, X = transformar(x)
        magnitud = np.abs(X)
        fase = np.angle(X)
        umbral = max(float(np.max(magnitud)) * 1e-3, 1e-12)
        fase_visible = np.where(magnitud >= umbral, fase, np.nan)

        ejes[fila, 0].plot(frecuencia, magnitud, lw=1.2)
        ejes[fila, 0].set_title(f"{nombre}: magnitud", loc="left", fontsize=10)
        ejes[fila, 0].set_ylabel("|X(f)| (amplitud·s)")
        ejes[fila, 0].grid(True, alpha=0.3)
        ejes[fila, 1].plot(frecuencia, fase_visible, lw=1.2)
        ejes[fila, 1].set_title(f"{nombre}: fase", loc="left", fontsize=10)
        ejes[fila, 1].set_ylabel("Fase (rad)")
        ejes[fila, 1].set_ylim(-np.pi - 0.15, np.pi + 0.15)
        ejes[fila, 1].set_yticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
        ejes[fila, 1].set_yticklabels(["−π", "−π/2", "0", "π/2", "π"])
        ejes[fila, 1].grid(True, alpha=0.3)

    for ax in ejes[-1, :]:
        ax.set_xlabel("Frecuencia (Hz)")
        ax.set_xlim(0, 100)
    fig.suptitle("Espectros de magnitud y fase (FFT de un solo lado)", fontsize=15)
    fig.tight_layout()
    fig.savefig(carpeta / "02_espectros_magnitud_fase.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)


def graficar_propiedades(pulso: np.ndarray, senoide: np.ndarray,
                         carpeta: Path) -> tuple[float, float]:
    a, b = 1.25, -0.60
    combinacion = a * pulso + b * senoide
    frecuencia, X_pulso = transformar(pulso)
    _, X_senoide = transformar(senoide)
    _, X_combinacion = transformar(combinacion)
    X_lineal = a * X_pulso + b * X_senoide
    error_lineal = float(
        np.max(np.abs(X_combinacion - X_lineal))
        / max(np.max(np.abs(X_combinacion)), np.finfo(float).eps)
    )

    tau = 0.125
    muestras = int(round(tau * FS))
    pulso_retrasado = np.roll(pulso, muestras)
    _, X_retrasado = transformar(pulso_retrasado)
    fase_medida = np.angle(X_retrasado * np.conj(X_pulso))
    fase_teorica = np.angle(np.exp(-2j * np.pi * frecuencia * tau))
    mascara = np.abs(X_pulso) > 0.02 * np.max(np.abs(X_pulso))
    error_desplazamiento = float(
        np.max(np.abs(np.angle(np.exp(1j * (fase_medida[mascara]
                                            - fase_teorica[mascara])))))
    )

    pulso_comprimido = ((np.abs(2.0 * T) <= 0.20)).astype(float)
    f_base, X_base = transformar(pulso)
    f_comprimido, X_comprimido = transformar(pulso_comprimido)

    fig, ejes = plt.subplots(1, 3, figsize=(16, 5.2))
    ventana = (T >= -0.35) & (T <= 0.35)
    ejes[0].plot(T[ventana], pulso[ventana], label="Pulso x(t)", alpha=0.75)
    ejes[0].plot(T[ventana], senoide[ventana], label="Senoide y(t)", alpha=0.75)
    ejes[0].plot(T[ventana], combinacion[ventana], label="1.25x(t) − 0.60y(t)", lw=1.4)
    ejes[0].set_title("Linealidad", fontsize=11)
    ejes[0].set_xlabel("Tiempo (s)")
    ejes[0].set_ylabel("Amplitud")
    ejes[0].legend(fontsize=8)
    ejes[0].grid(True, alpha=0.3)
    ejes[0].text(0.03, 0.03, f"Error relativo FFT: {error_lineal:.2e}",
                 transform=ejes[0].transAxes, fontsize=8,
                 bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})

    ejes[1].plot(frecuencia[mascara], fase_medida[mascara], ".",
                 ms=3, label="Diferencia medida")
    ejes[1].plot(frecuencia, fase_teorica, lw=1.1,
                 label="−2πfτ (envuelta)")
    ejes[1].set_xlim(0, 40)
    ejes[1].set_ylim(-np.pi - 0.2, np.pi + 0.2)
    ejes[1].set_title("Desplazamiento temporal (τ = 0.125 s)", fontsize=11)
    ejes[1].set_xlabel("Frecuencia (Hz)")
    ejes[1].set_ylabel("Diferencia de fase (rad)")
    ejes[1].legend(fontsize=8)
    ejes[1].grid(True, alpha=0.3)
    ejes[1].text(0.03, 0.03, f"Error máx.: {error_desplazamiento:.2e} rad",
                 transform=ejes[1].transAxes, fontsize=8,
                 bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})

    ejes[2].plot(f_base, np.abs(X_base), label="x(t): ancho ≈ 0.40 s")
    ejes[2].plot(f_comprimido, np.abs(X_comprimido),
                 label="x(2t): ancho ≈ 0.20 s")
    ejes[2].set_xlim(0, 30)
    ejes[2].set_title("Escalamiento en frecuencia", fontsize=11)
    ejes[2].set_xlabel("Frecuencia (Hz)")
    ejes[2].set_ylabel("|X(f)| (amplitud·s)")
    ejes[2].legend(fontsize=8)
    ejes[2].grid(True, alpha=0.3)

    fig.suptitle("Comprobación de propiedades de Fourier", fontsize=15)
    fig.tight_layout()
    fig.savefig(carpeta / "03_propiedades_fourier.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    return error_lineal, error_desplazamiento


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Genera señales, espectros y comprobaciones de Fourier."
    )
    parser.add_argument(
        "--salida", type=Path, default=Path(__file__).resolve().parent / "resultados",
        help="Carpeta donde se guardarán las gráficas (por defecto: ./resultados).",
    )
    args = parser.parse_args()
    args.salida.mkdir(parents=True, exist_ok=True)

    pulso, escalon, senoide = crear_senales(T)
    graficar_senales(pulso, escalon, senoide, args.salida)
    graficar_espectros([
        ("Pulso rectangular", pulso),
        ("Función escalón", escalon),
        ("Senoide de 50 Hz", senoide),
    ], args.salida)
    error_lineal, error_desplazamiento = graficar_propiedades(
        pulso, senoide, args.salida
    )

    print(f"Gráficas guardadas en: {args.salida.resolve()}")
    print(f"Error relativo de linealidad: {error_lineal:.3e}")
    print(f"Error máximo de fase en desplazamiento: {error_desplazamiento:.3e} rad")
    print("Nota: la FFT analiza muestras finitas; los bordes de la ventana afectan "
          "el espectro, especialmente el de la función escalón.")


if __name__ == "__main__":
    main()
