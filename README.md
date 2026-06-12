# Sistema de Saturno — Simulação Gravitacional N-Corpos

Simulação interativa do sistema de Saturno que roda inteiramente no navegador, em um único arquivo HTML. Todos os 10 corpos (Saturno + 9 luas) se atraem mutuamente pela **Lei da Gravitação Universal de Newton**, e as órbitas evoluem em tempo real a partir de dados astronômicos reais (JPL/NASA).

**▶ Acesse a simulação:** https://guilhermedrummondcouto.github.io/saturno/

> Substitua `SEU_USUARIO` acima pelo seu nome de usuário do GitHub depois de publicar.

## Destaques

- **Corpos simulados:** Saturno e suas 9 principais luas — Mimas, Encélado, Tétis, Dione, Reia, Titã, Hipérion, Jápeto e Febe (esta última em órbita retrógrada real).
- **Física N-corpos completa:** força gravitacional calculada entre todos os pares de corpos, com aceleração, velocidade e posição atualizadas continuamente.
- **Distâncias em escala real:** 1 unidade de cena = 1.000 km.
- **Anéis com estrutura real:** anéis D, C, B, divisão de Cassini, anel A, falhas de Encke e Keeler e anel F, com sombra do planeta projetada sobre eles.
- **Interação:** controle da velocidade do tempo (1 min/s a ~50 dias/s), pausa, trilhas orbitais, órbitas de referência, zoom/órbita/panorâmica da câmera, seleção de qualquer corpo e modo de acompanhamento (a câmera segue a lua escolhida).
- **Detalhes físicos extras:** acoplamento de maré das luas (sempre a mesma face voltada para Saturno), tombamento caótico de Hipérion, achatamento real de Saturno (0,098) e rotação de 10,7 h.

## Física utilizada

A aceleração de cada corpo *i* é a soma das atrações de todos os outros corpos *j*:

```
aᵢ = Σⱼ G·mⱼ·(rⱼ − rᵢ) / |rⱼ − rᵢ|³        (j ≠ i)
```

com G = 6,6743 × 10⁻²⁰ km³·kg⁻¹·s⁻², massas em kg e distâncias em km.

- **Integrador:** leapfrog simplético (*kick–drift–kick*) com passo Δt = 60 s. Integradores simpléticos não acumulam erro sistemático de energia, o que mantém as órbitas estáveis em longos períodos simulados.
- **Condições iniciais:** elementos orbitais reais de cada lua (semieixo maior, excentricidade, inclinação, nodo, periastro e anomalia) convertidos em vetores de posição e velocidade, no referencial do baricentro do sistema (momento total nulo).
- **Dados:** massas, raios e elementos orbitais aproximados do JPL/NASA Solar System Dynamics.

## Resultados da validação

A física foi validada por um script independente (`teste-validacao.js`, executável com Node.js: `node teste-validacao.js`) que replica exatamente o integrador da página e simula **600 dias (864.000 passos)**, medindo períodos orbitais, conservação de energia e estabilidade dos semieixos:

| Lua      | P Kepler (d) | P medido (d) | Erro %  | Δa %   | Status |
|----------|--------------|--------------|---------|--------|--------|
| Mimas    | 0,944        | 0,944        | 0,0007  | 0,0067 | OK     |
| Encélado | 1,371        | 1,370        | 0,0009  | 0,0047 | OK     |
| Tétis    | 1,888        | 1,888        | 0,0010  | 0,0054 | OK     |
| Dione    | 2,738        | 2,737        | 0,0002  | 0,0076 | OK     |
| Reia     | 4,519        | 4,518        | 0,0119  | 0,0157 | OK     |
| Titã     | 15,946       | 15,945       | 0,0032  | 0,0059 | OK     |
| Hipérion | 21,281       | 21,249       | 0,1523  | 1,9212 | OK     |
| Jápeto   | 79,338       | 79,199       | 0,1753  | 0,2032 | OK     |
| Febe     | 550,382      | 548,678      | 0,3096  | 0,3595 | OK     |

**Deriva máxima da energia total: 2,5 × 10⁻⁹ (relativa)** ao longo dos 600 dias — confirmando a estabilidade do integrador simplético.

A variação de ~1,9% no semieixo de Hipérion não é erro numérico: é física real, causada pela ressonância orbital 4:3 com Titã (Hipérion tem dinâmica caótica conhecida).

## Controles

| Ação                     | Como fazer                                   |
|--------------------------|----------------------------------------------|
| Orbitar a câmera         | Arrastar com o botão esquerdo                |
| Zoom                     | Roda do mouse ou gesto de pinça              |
| Mover (panorâmica)       | Arrastar com o botão direito (ou Shift)      |
| Selecionar um corpo      | Clique no corpo, no rótulo ou na lista       |
| Seguir um corpo          | Duplo clique, botão ◉ na lista ou "Seguir"   |
| Pausar / continuar       | Barra de espaço ou botão ⏸                   |
| Velocidade da simulação  | Controle deslizante no painel                |
| Visão geral              | Botão "Visão geral" ou tecla R               |

## Como executar localmente

Basta abrir o `index.html` em qualquer navegador moderno (é necessária conexão com a internet apenas para baixar a biblioteca Three.js do CDN).

## Créditos e licença

- Renderização 3D: [Three.js](https://threejs.org) (licença MIT).
- Dados orbitais: NASA/JPL Solar System Dynamics (valores aproximados).
- Texturas do planeta, dos anéis e das luas geradas proceduralmente em canvas.
- Código criado com [Claude](https://claude.com) (Anthropic).

Este projeto é distribuído sob a licença MIT — sinta-se livre para usar, modificar e compartilhar.
