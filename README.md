# Fly-by
Repository in portuguese.

##  Instalação
Apesar do código aparentar usar o CMake, isso é só por eu ter utilizado o CLion no desenvolvimento, e ser mais fácil configurar
a IDE com o CMake. Na prática, o código foi desenvolvido para ser compilado com o `gcc`. (A propósito, muito obrigado JetBrains por
fornecer uma versão gratuita do CLion para estudantes 🥹)

Para os gráficos, fui utilizado o [Julia](https://julialang.org/), que pode ser baixado no site oficial, ou usando:
```shell
curl -fsSL https://install.julialang.org | sh
```

Assim que o Julia estiver instalado, é possível adicionar os pacotes utilizados na construção dos gráficos a partir do documento
`install.jl`, rodando o comando:
```shell
julia install.jl
```

Com isso, o ambiente deve estar pronto para rodar a simulação e obter os resultados 😄

##  Simulação
A simulação do fly-by usando o problema dos dois corpos restrito é feita através do arquivo `fly_by_pr2c.c`. Para compilar ele você pode
rodar:
```shell
gcc fly_by_pr2c.c -lm -o fly_by_pr2c
```

É importante pontuar que esse código pode não funcionar no Windows, já que algumas bibliotecas utilizadas são específicas
de sistema Unix (como Linux e MacOS).

Para o problema de três corpos, o processo é semelhante, mas com o arquivo `fly_by_pr3c.c`:
```shell
gcc fly_by_pr3c.c -lm -o fly_by_pr3c
```

Tendo o executável compilado, é possível apenas executá-lo a fim de obter a lista de parâmetros de entrada. Algo como (não copie esse prompt):
```shell
./fly_by_pr2c
Use: ./fly_by_pr2c <test_name> <x_init_factor> <velocity_infinity> <b_min_factor> <b_max_factor> <max_time> <dt>
- <test_name>: Nome do teste, e da pasta onde a saída será salva.
- <x_init_factor>: Fator multiplicando R_Marte na definição de x(0). No trabalho usamos um valor igual a 20. Também é usado no critério de parada.
- <velocity_infinity>: Velocidade da sonda no infinito, em metros por segundo. No trabalho usamos o valor de 3000 m/s.
- <b_min_factor>: Fator mínimo usado na definição do intervalo de valores para o parâmetro de impacto. No trabalho usamos um valor igual a 2. Não recomendo tomar um valor menor do que 2.
- <b_max_factor>: Fator máximo usado na definição do intervalo de valores para o parâmetro de impacto. No trabalho usamos um valor igual a 6. Não recomendo tomar um valor maior do que 10.
- <max_time>: Critério de parada de emergência. É o tempo máximo que pode ser gasto com a integração antes dela ser abortada, sem segundos. No trabalho foi utilizado 1,8e5 segundos.
- <dt>: Passo temporal utilizado na integração, em segundos. Não deve ser muito grande já que é usado o método de Euler. No trabalho foi utilizado 0,001 segundos.
```

Para as simulações apresentadas no relatório da matéria de Simulação em Física, a execução foi feita com o comando:
```shell
 ./fly_by_pr2c simul 50 2600 -10 10 1e10 0.001
```

De maneria que os dados da simulação foram salvos na pasta `simul`. Uma execução bem sucedida vai resultar em uma saída como abaixo:
```shell
./fly_by_pr2c simul 50 2600 -10 10 1e10 0.001
Rodando o teste...
Condições iniciais definidas: 
         Raio de Marte utilizado: 3.3895e+06 metros 
         Massa de Marte utilizada: 6.4171e+23 kg 
         Valor de x(0): -1.6948e+08 metros
         Valor de y(0) pertencente ao intervalo [-3.3895e+07 m; 3.3895e+07 m], com passo igual a 2.8364e+05 metros
         Valor de vx(0): 2.6954e+03 metros por segundo
         Valor de vy(0): 0.0000e+00 metros por segundo
         Tempo máximo de integração: 1.0000e+10 segundos
         Passo de integração: 0.0010 s
Os dados serão salvos na pasta: 'simul'
Pasta do problema de 2 corpos: 'simul/pr2c'

Realizando simulações ... 
[##########################################################] 100.00%, Total time: 5 minutos e 20 segundos                                                                         

Salvando os dados globais em: 'simul/global.csv'
Simulação concluída =D

```

Após rodar a simulação de 2 corpos, rode a simulação do problema de 3 corpos. As configurações são semelhantes, mas é necessário definir também
a posição angular inicial de Marte. A execução e feita como a seguir:

```shell
./fly_by_pr3c simul 50 -0.01 2600 -10 10 1e10 0.001
```

A entrada aqui segue a mesma lógica, e pode ser facilmente verificada usando:
```shell
./fly_by_pr3c
Use: ./fly_by_pr3c <test_name> <x_init_factor> <velocity_infinity> <b_min_factor> <b_max_factor> <max_time> <dt>
- <test_name>: Nome do teste, e da pasta onde a saída será salva.
- <r_factor>: Fator multiplicando R_Marte que define o raio de influência do planeta. No trabalho usamos um valor igual a 50.
- <mars_init_angle>: Ângulo inicial de Marte no sistema de coordenadas cartesiano referenciado no Sol, em graus. No trabalho usamos um valor igual a -0.01
- <velocity_infinity>: Velocidade da sonda no infinito, em metros por segundo. No trabalho usamos o valor de 2600 m/s.
- <b_min_factor>: Fator mínimo usado na definição do intervalo de valores para o parâmetro de impacto. No trabalho usamos um valor igual a -10
- <b_max_factor>: Fator máximo usado na definição do intervalo de valores para o parâmetro de impacto. No trabalho usamos um valor igual a 10.
- <max_time>: Critério de parada de emergência. É o tempo máximo que pode ser gasto com a integração antes dela ser abortada, sem segundos. No trabalho foi utilizado 10e10 segundos.
- <dt>: Passo temporal utilizado na integração, em segundos. Não deve ser muito grande já que é usado o método de Euler. No trabalho foi utilizado 0,001 segundos.
```

É importante notar que o nome do teste deve ser igual; além disso, as configurações devem ser semelhantes para mante a proposta do trabalho. Uma
outra consideração é que a simulação de 3 corpos teve o código projetado para rodas apenas após a de 2 corpos, então, mudar a ordem da execução
pode resultar em erros na criação das pastas.

A saída esperada para essa simulação é algo como:
```shell
./fly_by_pr3c simul 50 -0.01 2600 -10 10 1e10 0.001
Rodando o teste...
Condições iniciais definidas: 
         Massa do Sol utilizada: 1.9885e+30 kg 
         Raio de Marte utilizado: 3.3895e+06 metros 
         Massa de Marte utilizada: 6.4171e+23 kg 
         Raio da órbita de Marte utilizada: 2.2794e+11 metros
         Valor raio de influência da esfera é de: 1.6948e+08 metros
         Valor do parâmetro de impacto pertencente ao intervalo [-3.3895e+07 m; 3.3895e+07 m], com passo igual a 2.8364e+05 metros
         Valor do módulo da velocidade inicial da sonda: 2.6954e+03 metros por segundo
         Tempo máximo de integração: 1.0000e+10 segundos
         Passo de integração: 0.0010 s
Pasta do problema de 3 corpos: 'simul/pr3c'

Realizando simulações ... 
[##########################################################] 100.00%, Total time: 20 minutos e 6 segundos                                                                         

Salvando os dados globais em: 'simul/global.csv'
Simulação concluída =D

```

## Gráficos
Tendo os dados da simulação, é possível obter os gráficos ao rodar o código
```shell
julia graphics.jl
```

Com os devidos pacotes instalados pelo script `install.jl`, o `graphics.jl` vai montar os gráficos da simulação, 
salvando eles na pasta `<test_name>/results/`. Abaixo tem um exemplo do output do script.
```shell
julia graphics.jl
Executando script gráficos em '/Users/gabrielferreira/Fly-by' ...
Por favor, selecione um desses seguintes projetos para iniciar o processamento:
        [1]      - Project: simul
Identificador do projeto: 1
Por favor, insira a velocidade no infinito (em m/s): 2600
Inicializando processamento dos dados globais do projeto 'simul' ...
Inicializando processamento das trajetórias do projeto 'simul' ...
Progress: 100%|████████████████████████████████████████████████████████████| Time: 0:00:21
```

Com isso, serão geradas as figuras apresentadas no arquivo `relatorio.pdf`, além de uma sequência de snapshots das trajetórias,
que são salvos em `<test_name>/results/snapshots`.