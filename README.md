# Simulador Didático de Pneumonia Bacteriana

Modelo matemático em Python que simula a evolução de uma infecção bacteriana no pulmão ao longo do tempo, com e sem antibiótico, e gera uma **imagem com os gráficos do resultado**.

> ⚠️ **Projeto educacional.** Os parâmetros são ilustrativos e o modelo é simplificado. Esta ferramenta **não** serve para diagnóstico, previsão de evolução clínica nem decisão de tratamento.

## 🎯 O que o projeto faz

A pneumonia bacteriana pode surgir como complicação de infecções respiratórias mal tratadas. Este projeto traduz essa ideia em um sistema de equações diferenciais (EDOs) e compara três cenários:

1. **Sem tratamento**
2. **Antibiótico com bactéria sensível**
3. **Antibiótico com bactéria resistente**

O resultado é salvo como imagem (`simulacao_pneumonia.png`) e como arquivo de dados (`resultados.json`).

## 📊 O que a imagem mostra

- **Gráfico superior:** carga bacteriana (log10 CFU/g de tecido pulmonar) ao longo de 96 horas, nos três cenários, com uma linha marcando o início do antibiótico.
- **Gráfico inferior:** nível de células imunes (neutrófilos/macrófagos), onde 1 é o valor basal.

## ⚙️ Como funciona

O modelo acompanha três quantidades:

| Variável | Significado |
|---|---|
| `B` | Carga bacteriana (CFU por grama de tecido) |
| `I` | Células imunes (1 = nível basal) |
| `A` | Concentração do antibiótico (múltiplos da CIM/MIC) |

- **Bactérias:** crescem de forma logística (limitadas por uma capacidade de suporte), são eliminadas pelas células imunes com fagocitose saturável e são mortas pelo antibiótico (efeito tipo Emax, que depende da concentração e é reduzido pela resistência).
- **Células imunes:** têm produção e morte basais e são recrutadas conforme a carga bacteriana, até um teto máximo.
- **Antibiótico:** soma de doses repetidas, cada uma com decaimento exponencial (meia-vida). O tratamento começa em 24 h, com doses a cada 8 h.
- **Resistência:** valor de 0 (sensível) a 1 (totalmente resistente). O cenário resistente usa 0,95.
- **Solução numérica:** `scipy.integrate.solve_ivp` (método LSODA).

## 🚀 Como usar

Requer Python 3.9 ou superior.

```bash
pip install -r requirements.txt
python simulador_pneumonia.py
```

Serão gerados:

- `simulacao_pneumonia.png`: imagem com os gráficos
- `resultados.json`: dados da simulação (tempo, log10 da carga bacteriana e células imunes de cada cenário)

Para testar outros cenários, altere os valores do dicionário `P` no início do arquivo (taxa de crescimento, dose, intervalo entre doses, início do tratamento etc.).

## 🚧 Limitações

- Modelo simplificado, que não representa a complexidade de um pulmão real nem a variação entre pacientes.
- Parâmetros ilustrativos, não calibrados com dados clínicos.
- Considera uma única população bacteriana e um único antibiótico.
- Não inclui comorbidades, dano tecidual, diferentes espécies de bactéria nem evolução da resistência durante o tratamento.
- Não substitui avaliação médica.

## 🗺️ Próximos passos

- Visualização 3D do pulmão (Three.js) lendo o `resultados.json`
- Substituir os parâmetros ilustrativos por valores de artigos científicos
- Incluir dano tecidual e cenários com atraso no início do tratamento
- Analisar a sensibilidade do modelo aos parâmetros

## 📚 Referências

Os parâmetros atuais são ilustrativos. Ao calibrá-los com a literatura, liste aqui os artigos utilizados.

## 🤖 Transparência

Este projeto foi desenvolvido com auxílio de IA (Claude, da Anthropic) e adaptado/estudado por mim como parte do meu aprendizado em bioinformática.

## 👩‍💻 Autora

Taciane · Biomedicina + Programação

[LinkedIn](COLE-SEU-LINK-AQUI) · [GitHub](COLE-SEU-LINK-AQUI)

## 📜 Licença e como citar

Distribuído sob a licença MIT (veja o arquivo `LICENSE`). Para citar, use o botão "Cite this repository" do GitHub (gerado a partir do `CITATION.cff`) ou o DOI do Zenodo.
