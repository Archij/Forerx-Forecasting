
//+------------------------------------------------------------------+
//|                                             Maximum Likeness.mq4 |
//|                  Copyright © 2014, by Артур Михайлович Степченко |
//|                                                                  |
//+------------------------------------------------------------------+
#property copyright "Copyright © 2014, by Артур Михайлович Степченко"
#property link      ""

#import "stdlib.ex4"
string ErrorDescription(int error_code);


extern string Indicators_=" Настройки индикатора";
/* Здесь прописываем настройки Вашего индикатора\индикаторов */
/* Стандартные переменные для шаблона НЕ ИЗМЕНЯТЬ!!! */
extern string trade_="Настройки торговли";
extern int Magic=666; // Магический номер
extern int OpenSlippage=1; // Проскальзывание
extern int CloseSlippage=3; // Проскальзывание 
extern int M = 24;
extern int P = 1;
extern int Risk = 1;

void Transpose(double &Matrix[][], double &transposeMatrix[][]) //Matrix transpose
{
  for (int i = 0; i<ArrayRange(Matrix, 0); i++)
  {
       for (int j = 0; j<ArrayRange(Matrix, 1); j++)
       {
            transposeMatrix[j][i]=Matrix[i][j];
       }
  }
}

void Multiply(double &Matrix1[][],double &Matrix2[][], double &Result[][]) //Matrix multiply
{
  int m1rows = ArrayRange(Matrix1, 0);
  int m1cols = ArrayRange(Matrix1, 1);
  int m2cols = ArrayRange(Matrix2, 1);
  
  for (int i=0; i<m1rows; i++)
  {
      for (int j=0; j<m2cols; j++)
      {
        for (int k=0; k<m1cols; k++)
        {
            Result[i][j] = Result[i][j] + Matrix1[i][k] * Matrix2[k][j];
        }
       }
  }
}

void Inverse(double &Matrix[][], double &invMatrix[][]) //Inverse matrix
{
    int i,j,k,n;
    double temp;
    
    n = ArrayRange(Matrix, 0); //matrix dimension
    
    for (i=0;i<n;i++) //creation of unit matrix
    {
        invMatrix[i][i]=1;
        for(j=0;j<n;j++)
        {
            if(i!=j)
            {
               invMatrix[i][j]=0;
            }
        }
    }
    for(i=0;i<n;i++)
    {
        temp=Matrix[i][i];
        if(temp==0)
        {
           break;
        }
        for(j=0;j<n;j++)
        {
            Matrix[i][j]=Matrix[i][j]/temp;
            invMatrix[i][j]=invMatrix[i][j]/temp;
        }
            for(k=0;k<n;k++)
            {
                if (k!=i)
                {
                    temp=Matrix[k][i];
                    for(j=0;j<n;j++)
                    {
                        Matrix[k][j]=Matrix[k][j]-temp*Matrix[i][j];
                        invMatrix[k][j]=invMatrix[k][j]-temp*invMatrix[i][j];
                    }
 
                }
            }
    }
  
}

int Extrapolation()
{
   int Sig = 0;
   
   int T = Bars-1;
  
   int N = M; 
  
  double XNM[][1]; //вектор длины M, лежащий внутри исходного X(t) началом которого является момент времени t=N
  
  ArrayResize(XNM, M);
  
  for (int t = 0; t<M; t++)
  {
      XNM[t][0] = Close[N-t];
  }
  
  double L[]; //Вектор значений моделей коэффициентов линейной корреляции (Вектор подобия)
  
  ArrayResize(L, T-N);
  
  double XiM[];
  
  ArrayResize(XiM, M);
  
  double XNMaverage = 0;
  double XiMaverage = 0;
  double cov = 0; //Ковариация исходных векторов
  //Дисперсии
  double DXNM = 0;
  double DXiM = 0;
  int J = 0;

  for (int i=Bars; i>N; i--)
  {
      
     if ((N >=1 && N <= T-1) && (M >=1 && M <= T-1) && (M+N)<T && (M+J)<T)
     {
     
        for (int iM = 0; iM<M; iM++)
        {
             XiM[iM] = Close[i-iM];
        }
     
        for (int xnmav = 0; xnmav<M; xnmav++)
        {
            XNMaverage = XNMaverage + XNM[xnmav][0];
        }
        
         XNMaverage = XNMaverage / M;
         
        for (int ximav = 0; ximav<M; ximav++)
        {
            XiMaverage = XiMaverage + XiM[ximav];
        }
        
        XiMaverage = XiMaverage / M;
        
        for (int c = 0; c<M; c++)
        {
            cov = cov + (XNM[c][0]-XNMaverage)*(XiM[c]-XiMaverage);
        }
        
        for (int d = 0; d<M; d++)
        {
            DXNM = DXNM + MathPow((XNM[d][0]-XNMaverage),2);
            DXiM = DXiM + MathPow((XiM[d]-XiMaverage),2);
        }

        L[J] = MathAbs(cov / MathSqrt(DXNM*DXiM));
        
        J = J+1;
 
     }

  }
  
  int Lindex = ArrayMaximum(L,WHOLE_ARRAY,0);
  
  double XimaxM[][2];
  
  ArrayResize(XimaxM, M);
  
  for (int imax = 0; imax<M; imax++)
  {
        XimaxM[imax][0] = Close[Bars-Lindex-imax];
        XimaxM[imax][1] = 1;
  }
  
  double XimaxMT[2][24];
  
  //Matrix transpose
  Transpose(XimaxM,XimaxMT);
 
  double result1[2][2];
  
  //Matrix multiply
  Multiply(XimaxMT, XimaxM, result1);

  double inverse[2][2];
   
  //Matrix inverse
  Inverse(result1, inverse);
  
  double result2[2][24];
  
  //Matrix multiply
  Multiply(inverse, XimaxMT, result2);
  
  double A[2][1];  //матрица линейных коэффициентов размерностью 2х1
  
  //Matrix multiply
  Multiply(result2, XNM, A);
  
  double XTP[1][1];
  
  double Xprim[2][1];
  
  ArrayResize(Xprim, P);
  
  for (int pr = 0; pr<P; pr++)
  {
        Xprim[0][pr] = Close[Bars-Lindex-M-pr];
        Xprim[1][pr] = 1;
  }
  
  double AT[1][2];
  
  Transpose(A,AT);
  
    //Matrix multiply
  Multiply(AT, Xprim, XTP);
  
  double Price = XTP[0][0];
  
  Alert(Price);
  
  if (Open[0] < Price)
  {
     Sig = 1;
  }
  else if (Open[0] > Price)
  {
     Sig = 2;
  }

   return(Sig);
}

int index = 1;

//+------------------------------------------------------------------+
//| expert initialization function                                   |
//+------------------------------------------------------------------+
int buy_ticket = -1;
int sell_ticket = -1;
int err = 0;
bool order_Close = false;
static int prevtime = 0;
double Lot = 0;

int init()
  {

     prevtime = iTime(Symbol(),0,0);

     return(0);
  }
//+------------------------------------------------------------------+
//| expert deinitialization function                                 |
//+------------------------------------------------------------------+
int deinit()
  {

   return(0);
  }
//+--------------------------------------------------- ---------------+
//| expert start function                                            |
//+------------------------------------------------------------------+
int start()
  {

    if (iTime(Symbol(), 0, 0) == prevtime)
    {
       return(0);
    }

    prevtime = iTime(Symbol(),0,0); 

    RefreshRates();
    

    if (buy_ticket != -1)
    {
            
              order_Close = OrderClose(buy_ticket,Lot,Bid,CloseSlippage,Green);
           
              if (order_Close == false)
              {
                  err = GetLastError();
                  Alert("Buy order closing failed with error: ",ErrorDescription(err));
                  err = 0;
              }
              else
              {
                   buy_ticket = -1; 
              }
                        

      }

      if (sell_ticket != -1)
      {
            
                  order_Close = OrderClose(sell_ticket,Lot,Ask,CloseSlippage,Red);
          
                  if (order_Close == false)
                  {
                    err = GetLastError();
                    Alert("Sell order closing failed with error: ",ErrorDescription(err));
                    err = 0;
                  }
                  else
                  {
                    sell_ticket = -1;
                  }
              
                      
    }

   

    if (Extrapolation()==1 && buy_ticket == -1)
    {  
               Lot = (Risk*AccountFreeMargin())/100 * AccountLeverage() * Close[0] / 100000;
           
               if (Lot > MarketInfo(Symbol(), MODE_MAXLOT))
               {
                  Lot = MarketInfo(Symbol(), MODE_MAXLOT);
               }
               else if (Lot < MarketInfo(Symbol(), MODE_MINLOT))
               {
                  Lot = MarketInfo(Symbol(), MODE_MINLOT);
               }
            
               buy_ticket = OrderSend(Symbol(), OP_BUY, Lot, Ask, OpenSlippage, 0, 0, "Buy Order", Magic, 0, Green);
               
               if (buy_ticket<0)
               {
                   err = GetLastError();
                   Alert("Buy order opening failed with error: ",ErrorDescription(err));
                   err = 0;
                   buy_ticket = -1;
               }   
           
             
    }
    
    if (Extrapolation()==2 && sell_ticket == -1)
    {
               Lot = (Risk*AccountFreeMargin())/100 * AccountLeverage() * Close[0] / 100000;
           
               if (Lot > MarketInfo(Symbol(), MODE_MAXLOT))
               {
                  Lot = MarketInfo(Symbol(), MODE_MAXLOT);
               }
               else if (Lot < MarketInfo(Symbol(), MODE_MINLOT))
               {
                  Lot = MarketInfo(Symbol(), MODE_MINLOT);
               }
               
               sell_ticket = OrderSend(Symbol(), OP_SELL, Lot, Bid, OpenSlippage, 0, 0, "Sell Order", Magic, 0, Red);
               
               if (sell_ticket<0)
               {
                  err = GetLastError();
                  Alert("Sell order opening failed with error: ",ErrorDescription(err));
                  err = 0;
                  sell_ticket = -1;

               }
        
    
    }
    
    return(0);
 }