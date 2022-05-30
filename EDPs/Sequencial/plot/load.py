import numpy as np
import pandas as pd
from typing import Iterable

# Funções que podem ser importadas ---------------------------------------------

# Exemplo de import:
# from load import read_np_from_csv, write_np_to_csv

# Exemplo de uso pode ser encontrado a baixo na função 'test()'

# Ler vetores do CSV
# A função causará um erro se o número de variáveis e o número de colunas for diferente

# Se você quiser ler só algumas colunas e descartar o resto, pode ser feito:
# >>> x, y, *outras_colunas_que_nao_me_importam = read_np_from_csv('test.csv')
# OU
# >>> x, *resto, valor = read_np_from_csv('test.csv')
# OU
# >>> x, y, *resto, h, p = read_np_from_csv('test.csv')
# E quaisquer outras combinações...
# O importante é ter no mínimo o mesmo número de colunas que o número de variáveis sem o asterisco

#x_, y_, valor_ = read_np_from_csv('test.csv')

# Garantir que funciona
#assert np.isclose(x, x_).all()
#assert np.isclose(y, y_).all()
#assert np.isclose(valor, valor_).all()

#print("Os arrays são iguais! As funções funcionam!")

# Escrever vetores no CSV
# Parâmetros: nome do arquivo, *vários arrays de numpy....., OPCIONALMENTE: columns=lista dos nomes

# Note que o ´columns=´ é necessário para especificar os nomes por conta dos inúmeros parâmetros posicionais.
# Os nomes podem ser omitidos. A chamada abaixo poderia ser somente:
# >>> write_np_to_csv('test.csv', x, y, valor)
# A única diferença é que o arquivo final teria as colunas como '0', '1', '2',
# mas não afetaria a leitura pela função ´read_np_from_csv´

#write_np_to_csv('test.csv', x, y, valor, columns=('x', 'y', 'valor'))

def split_strip(text, sep):
    """ Separar text por sep e remover espaços em volta """
    return [s.strip() for s in text.split(sep)]

def read_np_from_csv(filename: str) -> 'tuple[np.array, ...]':
    return pd.read_csv(filename).to_numpy().T

def write_np_to_csv(filename: str, *data, columns: 'Iterable[str]'=None) -> None:
    data_matrix = np.matrix(data).T
    df = pd.DataFrame(data_matrix, columns=columns)
    df.to_csv(filename, index=False)

# TESTES -----------------------------------------------------------------------

def test():
    from random import random

    print("Rodando testes...")

    # Valores aleatórios
    #x = np.random.random((10,))
    #y = np.random.random((10,))
    #valor = np.random.random((10,))

    x = read_np_from_csv("x.csv")
    y = read_np_from_csv("y.csv")
    z = read_np_from_csv("z.csv")
    t = read_np_from_csv("t.csv")
    bacteria = read_np_from_csv("bacteria1d.csv")    
    print(x)
    print('-------------------------')
    #print(y)
    print('-------------------------')
    #print(z)
    print('-------------------------')
    #print(bacteria)
    #print('1º Conjunto de valores em x:\n' + str(bacteria[:][0]))
    values = bacteria[:][1] 
    print(values)   
    print(np.size(x))
    r = np.hsplit(values, np.size(x))
    print('r[0] ' + str(r[0]))
    print('r ' + str(r))
    
    i = 0
    for v in r: 
        print('t=' + str(t[0][i]) + ':\n')
        #print(v)
        i+=1
    #print('2ª Linha:\n' + bacteria[1])


if __name__ == '__main__':
    test()