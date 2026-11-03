"""
Simulador didático de pneumonia bacteriana (modelo de EDOs)

Variáveis:
    B = carga bacteriana (CFU por grama de tecido pulmonar)
    I = nível de células imunes (neutrófilos/macrófagos), 1 = basal
    A = concentração do antibiótico (em múltiplos da CIM/MIC)

ATENÇÃO: parâmetros ilustrativos. Para algo sério, substitua por valores
de artigos científicos. Uso didático, NÃO clínico.
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ---------------- Parâmetros ----------------
P = dict(
    r=0.5,            # crescimento bacteriano (1/h)
    K=1e9,            # capacidade de suporte (CFU/g)
    kI=0.05,          # morte bacteriana por célula imune (1/h por unidade de I)
    h=1e7,            # saturação da fagocitose (CFU/g)
    s=0.1, d=0.1,     # produção e morte basal de células imunes (basal I = s/d = 1)
    aI=0.4,           # recrutamento imune estimulado pela bactéria (1/h)
    Imax=20.0,        # teto de células imunes
    kA=1.5,           # morte máxima pelo antibiótico (1/h)
    EC50=1.0,         # concentração com metade do efeito (x CIM)
    dose_peak=4.0,    # pico após cada dose (x CIM)
    half_life=2.0,    # meia-vida do antibiótico (h)
    dose_every=8.0,   # intervalo entre doses (h)
    t_start=24.0,     # início do tratamento (h)
    n_doses=9,
)


def antibiotic(t, p):
    """Concentração do antibiótico no tempo t (soma de doses com decaimento exponencial)."""
    if t < p["t_start"]:
        return 0.0
    ke = np.log(2) / p["half_life"]
    n = min(int((t - p["t_start"]) // p["dose_every"]) + 1, p["n_doses"])
    return sum(p["dose_peak"] * np.exp(-ke * (t - (p["t_start"] + k * p["dose_every"])))
               for k in range(n))


def rhs(t, y, p, treat, resistance):
    B, I = y
    B = max(B, 0.0)
    A = antibiotic(t, p) if treat else 0.0
    sat = B / (B + p["h"])
    # (1 - resistance): 0 = bactéria sensível, 1 = totalmente resistente
    dB = (p["r"] * B * (1 - B / p["K"])
          - p["kI"] * I * B / (1 + B / p["h"])
          - p["kA"] * (A / (A + p["EC50"])) * B * (1 - resistance))
    dI = p["s"] - p["d"] * I + p["aI"] * I * sat * (1 - I / p["Imax"])
    return [dB, dI]


def simulate(treat, resistance, p=P, t_end=96.0, B0=1e3):
    t_eval = np.linspace(0, t_end, 400)
    sol = solve_ivp(rhs, (0, t_end), [B0, 1.0], t_eval=t_eval,
                    args=(p, treat, resistance), method="LSODA", rtol=1e-8, atol=1e-8)
    return sol.t, np.maximum(sol.y[0], 1.0), sol.y[1]


if __name__ == "__main__":
    scenarios = {
        "Sem tratamento": (False, 0.0),
        "Antibiótico (bactéria sensível)": (True, 0.0),
        "Antibiótico (bactéria resistente)": (True, 0.95),
    }

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    export = {}
    for nome, (treat, res) in scenarios.items():
        t, B, I = simulate(treat, res)
        ax1.plot(t, np.log10(B), label=nome, lw=2)
        ax2.plot(t, I, label=nome, lw=2)
        export[nome] = {"t": t.tolist(), "log10_B": np.log10(B).tolist(), "I": I.tolist()}

    ax1.axvline(P["t_start"], color="gray", ls="--", alpha=.6)
    ax1.text(P["t_start"] + 1, 0.3, "início do antibiótico", color="gray")
    ax1.set_ylabel("Carga bacteriana (log10 CFU/g)")
    ax1.set_title("Simulação didática de pneumonia bacteriana")
    ax1.legend()
    ax1.grid(alpha=.3)
    ax2.set_ylabel("Células imunes (1 = basal)")
    ax2.set_xlabel("Tempo (horas)")
    ax2.grid(alpha=.3)
    plt.tight_layout()
    plt.savefig("simulacao_pneumonia.png", dpi=150)
    plt.show()

    # Esse JSON será a ponte para a visualização 3D (Three.js) depois
    with open("resultados.json", "w", encoding="utf-8") as f:
        json.dump(export, f)
    print("Gerados: simulacao_pneumonia.png e resultados.json")

   