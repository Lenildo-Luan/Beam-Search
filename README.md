# Scheduling a single machine with release dates and sequence dependent setup times
 
Para entender o problema, imagine uma fábrica que possui uma máquina que será utilziada para a criação de alguns produtos diferentes. Imagine que para a fabricação de cada um desses produtos seja necessário configurar a máquina, esta configuração pode levar um tempo diferente de acordo com o produto que foi fabricado anteriormente. Além disso, cada produto demora um determinado tempo para ser criado. Outro ponto importante, é que alguns produtos só estão disponíveis para serem criados a partir de um determinado tempo (ex: uma hora após o início do processo).

O objetivo deste problema é encontrar a alocação de recursos de forma a minimizar o a soma total do tempo de espera até o termino da fabricação do último produto.
 
## Formas de resolução
 
Existem dois métodos que podem ser seguidos para resolver problemas como o descrito anteriormente, os exatos ou os heurísticos.

Métodos exatos visam iterar sobre todo o conjunto de soluções de um problema para obter sua resposta. Por conta disso, eles retornam a solução ótima do problema, ou seja, a melhor dentre todas as soluções. Por outro lado, para problemas da classe NP, não é possível encontrar estas soluções em tempo computacional aceitável a partir de certo ponto. Para ilustrar, uma instância deste problema com 60 vértices possui uma quantidade de soluções possíveis semelhante a quantidade de átomos no universo!

Por outro lado, existem os métodos heurísticos. Eles tentam resolver estes problemas com um tempo computacional aceitável, porém sem garantir sua otimalidade. Assim como as meta-heurísticos, que tem a mesma finalidade, porém com o intuito de resolver vários tipos de problemas diferentes com um mesmo algoritmo. Por mais que a otimalidade não seja garantida, muitas heurísticas e meta-heurísticas conseguem encontrar resultados ótimos ou muito próximo deles para uma grande gama de problemas.

A implementação presente neste repositório busca resolver o problema de forma heurística através de um algorítimo denominado [Beam Search](https://www.sciencedirect.com/science/article/abs/pii/S0305054816300776?via%3Dihub) adaptado para resolver o presente problema.
 
## Instância do problema
 
As instâncias do problema encontram-se na pasta `./instances` e foram obtidas a partir da literatura existente. Você pode saber mais sobre elas [aqui](https://www.sciencedirect.com/science/article/abs/pii/S0305054816300776?via%3Dihub).
 
## Compilando e executando
 
Para compilar o programa, basta utilizar o comando `make` pelo terminal na pasta raiz, gerando o arquivo executável `bs`. Basta executar o arquivo gerado passando como parâmetro a instância a ser utilizada para ver o resultado.
 
``` bash
./bs instances/in02_001.dat b
```
 
## Resultados
 
## TODO
 
* Criar tabela com os resultados computacionais
