# Forerx Forecasting

Набор разных методов машинного обучения для прогнозирования нелинейных временных рядов из бирж Forex, а также разные технические индикаторы для Forex.

MATLAB функции

WRVM.m - взвешенный метод релевантных векторов.

PHI.m - имплементация ядровых преобразований, которые можно использовать для метода релевантных векторов или метода опорных векторов.

trainLRNN.m - подбор подходящих параметров для модели Layer Recurrent Neural Network (LRNN).

knn.m - подбор подходящих параметров для модели к-ближайших соседей (kNN).

Prepare_Train_Data.m - подготовка набора данных для входа моделей машинного обучения.

DifferentialPhaseSpace.m - реконструкция дифференциального фазового пространства для хаотичных временных рядов. Код написан по алгоритму из статьи: https://www.sciencedirect.com/science/article/pii/S0307904X07003526

HMA.m - скользящая средняя Хала.

JMA.m - скользящая средняя Юрика.

LWMA.m - линейно-взвешенная скользящая средняя.

NonLagMA.m - скользящая средняя "без лагов".

ZeroLagEMA.m -  экспоненциальная скользящая средняя с "нулевым лагом".

ZeroLagJMA.m - скользящая средняя Юрика c "нулевым лагом".

ZigZag.m - зиг-заг индикатор.

MQL4 функции

Maximum Likeness.mq4 - торговый робот, основанный на максимуме подобия.

RVM.mq4 - торговый робот, основанный на методе релевантных векторов.

Artificial Neural Networks.mq4 - торговый робот, основанный на нейронной сетью с прямой связью.
