
//+------------------------------------------------------------------+
//|                                   Artificial Neural Networks.mq4 |
//|                  Copyright © 2014, by Артур Михайлович Степченко |
//|                                                                  |
//+------------------------------------------------------------------+
#property copyright "Copyright © 2014, by Артур Михайлович Степченко"
#property link      ""

#import "stdlib.ex4"
string ErrorDescription(int error_code);


extern int Magic = 666; //Магический номер
extern bool AutoLot = 0; //0 - константный лот, 1 - динамичный лот
extern double Lot = 0.1; //Константный лот
extern int Risk = 10; //Риск в процентах от баланса при динамичном лоте
extern int TakeProfit = 0; //Тейк профит в пунктах, 0 - не используется
extern int StopLoss = 0; //Стоп лосс в пунктах, 0 - не используется
extern int OpenSlippage=3; // Проскальзывание цены открытия ордера
extern int CloseSlippage=3; // Проскальзывание цены закрытия ордера
extern int Transformation = 0; //Трансформация: 0 - лог, 1 - проц, 2 - разн
extern int Normalization = 0; //Нормализация: 0 - [-1;1], 1 - [0;1]
extern int Input_Nodes = 3; //Окно (глубина погружения)
extern double nu = 0.3; //Корректирующий коэффициент nu
extern double lamda = 0.2; //Корректирующий коэффициент lamda
extern double Emin = 0.05; //Минимальная ошибка нейроной сети

double W[]; //Вектор весов выходного слоя
double V[]; //Вектор весов входного слоя
double Wbias;
double Vbias[];
double close[];
double Mean;
double Standart_Deviation;
double x_min;
double x_max;
double O;
int Hidden_Nodes = 3;//Input_Nodes*2+1; //Число нейронов скрытого слоя


//+------------------------------------------------------------------+
//| expert initialization function                                   |
//+------------------------------------------------------------------+
int buy_ticket = -1;
int sell_ticket = -1;
int A = 0;
bool order_Close = false;
static int prevtime = 0;
int index = 1;
double Price[];


void Standardization()
{
    
    Mean = 0;
    
    for (int i=0;i<ArraySize(Price);i++)
    {
         Mean = Mean + Price[i];

    }
    
    Mean = Mean/ArraySize(Price);
    
    for (i=0;i<ArraySize(Price);i++)
    {
         Price[i] = Price[i] - Mean;
    }
     
    double Variance = 0;
    
    for (i=0;i<ArraySize(Price);i++)
    {
         Variance = Variance + MathPow(Price[i],2);     
    }
     
    Variance = Variance/ArraySize(Price);
    
    Standart_Deviation = MathSqrt(Variance);
    
    for (i=0;i<ArraySize(Price);i++)
    {
        Price[i] = Price[i] / Standart_Deviation;

    }
     
}

void Logarithm()
{

    int ind = Bars-3;
    for (int i = 0; i<=ind; i++)
    {
        Price[ind-i] = log(Close[i+1]/Close[i+2]);
    }
    
}

void Percentage()
{

    int ind = Bars-3;
    for (int i = 0; i<=ind; i++)
    {
        Price[ind-i] = Close[i+1]/Close[i+2];
    }
}

void Difference()
{

    int ind = Bars-3;
    for (int i = 0; i<=ind; i++)
    {
        Price[ind-i] = Close[i+1] - Close[i+2];
    }
}

void Normalization(int a, int b, double xmin, double xmax)
{


        for (int i=0;i<ArraySize(Price);i++)
        {
            Price[i]=(((Price[i]-xmin)*(b-a))/(xmax-xmin))+a;

   
        }

}


int ANN()
{
      int Sig = 0;
      
      /*
      double lim1 = -0.5;
      double lim2 = 0.5;
      
      for (int i = 0; i<Hidden_Nodes; i++)
      {
          W[i] = NormalizeDouble(lim1 + (lim2-lim1)*MathRand()/32767.0,2); 
      }
      
      for (i = 0; i<Input_Nodes*Hidden_Nodes; i++)
      {
          V[i] = NormalizeDouble(lim1 + (lim2-lim1)*MathRand()/32767.0,2); 
      }
      
      Wbias = NormalizeDouble(lim1 + (lim2-lim1)*MathRand()/32767.0,2);
      
      for (i = 0; i<Hidden_Nodes; i++)
      {
          Vbias[i] = NormalizeDouble(lim1 + (lim2-lim1)*MathRand()/32767.0,2); 
      }
      
      */
      
      W[0] = 0.25;
      W[1] = 0.50;
      W[2] = 0.75;
      
      Wbias = 1.00;
      
      V[0] = 0.5;
      V[1] = 0.6;
      V[2] = 0.2;
      
      V[3] = 1.3;
      V[4] = 0.6;
      V[5] = 0.1;
      
      V[6] = 1.5;
      V[7] = 0.5;
      V[8] = 1.0;
      
      Vbias[0] = 1.2;
      Vbias[1] = 0.5;
      Vbias[2] = 0.7;

      double E = Emin+1;
      
      double x[];
      ArrayResize(x,Input_Nodes,0);
      
      double y[];
      ArrayResize(y,Hidden_Nodes,0);
         
      double net[];
      ArrayResize(net,Hidden_Nodes,0);
      
      ArrayInitialize(net,0);
      
      int Epochs = ArraySize(Price) - Input_Nodes;
      double d;
      int epoch = 1;
      
      double nety = 0.0;
      int nindex = 0;
      int vindex = 0;
      int yindex = 0;

      double b0;
      double by[];
      ArrayResize(by,Hidden_Nodes,0);
      

      
      while (E>=Emin || epoch > 2000)
      {

          E = 0.0;
          
          for (int step = 0; step<Epochs; step++)
          {

               for (int xi = 0; xi<Input_Nodes; xi++)
               {
    
                   x[xi] = Price[xi+step];
               }
               
            
               d = Price[Input_Nodes+step];

               for (int ni = 0; ni<Hidden_Nodes; ni++)
               {
                  
                  for (int nj = 0; nj<Input_Nodes; nj++)
                  {
                      net[ni] = net[ni] + V[nj+nindex]*x[nj];

                  }
                  
                  net[ni] = net[ni] + Vbias[ni];
                  
                  nindex = nindex + Input_Nodes;
                     
               }
 
               nindex = 0;
               
               //Computing output vector y[y1, y2, ..., yWindow] on hidden layer 

               for (int yi = 0; yi<Hidden_Nodes; yi++)
               {
                    y[yi] = 2 / (1 + exp(-lamda*net[yi]))-1;
                    
               }

               ArrayInitialize(net,0);

               for (int ny = 0; ny<Hidden_Nodes; ny++)
               {
                   nety = nety + W[ny]*y[ny];
               }
               
               nety = nety + Wbias;

               //Computing the output vector O on output layer 
               //O = nety;
               
               O = 2 / (1 + exp(-lamda*nety))-1;

               nety = 0.0;

               E = E + 0.5 * MathPow(d-O,2);
               
               //Stepest gradient Descent
               
               //Error for output layer
               b0 = (d-O) * (lamda/2)  * (1 - O) * (1 + O);
               
               //b0 = (d-O);

               //Errors for hidden layers
               for (int bi = 0; bi < Hidden_Nodes; bi++)
               {
                   by[bi] = (lamda/2) * (1-y[bi]) * (1+y[bi]) *  W[bi] * b0;
                    
               }
           
               for (int wi = 0; wi < Hidden_Nodes; wi++)
               {
                 
                      W[wi] = W[wi] - nu * b0 * y[wi];

               }
               
               Wbias = Wbias - nu*b0;
               
               for (int vi = 0; vi<Hidden_Nodes; vi++)
               {
                  
                   for (int vj = 0; vj<Input_Nodes; vj++)
                   {

                       V[vj+yindex] = V[vj+yindex] - nu*by[vi]*x[vj+yindex];

                        
                   }
                   
                   Vbias[vi] = Vbias[vi] - nu*by[vi];
                   
                   yindex = yindex - Input_Nodes;
  
               }

               yindex = 0;
               

               
               
          }
          
          if (E >= Emin)
          {
             epoch = epoch + 1;
          }
      
      
      }
      
      //экстраполяция
      
      for (xi = 0; xi<Input_Nodes; xi++)
      {
               
          x[xi] = Price[xi+Epochs];
      }
               
       for (ni = 0; ni<Hidden_Nodes; ni++)
       {
                  
            for (nj = 0; nj<Input_Nodes; nj++)
            {
                net[ni] = net[ni] + V[nj+nindex]*x[nj];

            }
                  
            net[ni] = net[ni] + Vbias[ni];
                  
            nindex = nindex + Input_Nodes;
                     
       }
 
       nindex = 0;
               
       //Computing output vector y[y1, y2, ..., yWindow] on hidden layer 

       for (yi = 0; yi<Hidden_Nodes; yi++)
       {
            y[yi] = 2 / (1 + exp(-lamda*net[yi]))-1;
                    
       }

       ArrayInitialize(net,0);

       for (ny = 0; ny<Hidden_Nodes; ny++)
       {
           nety = nety + W[ny]*y[ny];
       }
               
       nety = nety + Wbias;

       //Computing the output vector O on output layer 
       //O = nety;
               
       O = 2 / (1 + exp(-lamda*nety))-1;

       nety = 0.0;
      
      //Normalize back

      //O=(((O-Price[ArrayMinimum(Price,WHOLE_ARRAY,0)])*(x_max-x_min))/(Price[ArrayMaximum(Price,WHOLE_ARRAY,0)]-Price[ArrayMinimum(Price,WHOLE_ARRAY,0)]))+x_min;
      
      //O = (O * Standart_Deviation)+Mean; //Standartize back
      
      Alert(O);
      Alert(epoch);
      
      if (Transformation == 0)
      {
          //O = exp(O);
         // O = O * Close[1];
      }
      if (Transformation == 1)
      {
          O = O * Close[1];
      }
      if (Transformation == 2)
      {
          O = O + Close[1];
      }
      
      if (O > Close[1])
      {
         Sig = 1;
      }
      
      else if (O < Close[1])
      {
         Sig = 2;
      }
      
      return(Sig);
}

int init()
  {

     prevtime = iTime(Symbol(),0,0);

     MathSrand(GetTickCount());
     
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
    


    ArrayResize(Price, Bars-2, 0);
    
    //Трансформация
   
    
    if (Transformation == 0)
    {
       //Logarithm();
    }
    else if (Transformation == 1)
    {
       Percentage();
    }
    else if (Transformation == 2)
    {
       Difference();
    }
    
    //Standardization();
    
    //Нормализация
    
    x_min=Price[ArrayMinimum(Price,WHOLE_ARRAY,0)];
    x_max=Price[ArrayMaximum(Price,WHOLE_ARRAY,0)];
    
    if (Normalization == 0)
    {
       //Normalization(-1, 1, x_min, x_max);
    }
    else if (Normalization == 1)
    {
       Normalization(0, 1, x_min, x_max);
    }
    
    ArrayResize(W,Hidden_Nodes,0);
    ArrayResize(V,Input_Nodes*Hidden_Nodes,0);
    ArrayResize(Vbias,Hidden_Nodes,0);
    
    ArrayResize(Price,20,0);
    
    Price[0] = 1;
    Price[1] = 0.94;
    Price[2] = 0.92;
    Price[3] = 0.83;
    Price[4] = 0.82;
    Price[5] = 0.87;
    Price[6] = 0.88;
    Price[7] = 0.86;
    Price[8] = 0.84;
    Price[9] = 0.82;
    Price[10] = 0.76;
    Price[11] = 0.79;
    Price[12] = 0.84;
    Price[13] = 0.88;
    Price[14] = 0.86;
    Price[15] = 0.88;
    Price[16] = 0.89;
    Price[17] = 0.82;
    Price[18] = 0.72;
    Price[19] = 0.68;
    
    
    
    
    if (buy_ticket != -1)
    {
            
              order_Close = OrderClose(buy_ticket,Lot,Bid,CloseSlippage,Green);
           
              if (order_Close == false)
              {
                  A = GetLastError();
                  Alert("Buy order closing failed with error: ",ErrorDescription(A));
                  A = 0;
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
                    A = GetLastError();
                    Alert("Sell order closing failed with error: ",ErrorDescription(A));
                    A = 0;
                  }
                  else
                  {
                    sell_ticket = -1;
                  }
              
                      
       }

/*
    if (OrderSelect(buy_ticket,SELECT_BY_TICKET)==true && OrderTakeProfit() != 0 && OrderClosePrice() >= OrderTakeProfit())
    {
         buy_ticket = -1;
    } 
    
    if (OrderSelect(sell_ticket,SELECT_BY_TICKET)==true && OrderTakeProfit() != 0 && OrderClosePrice() <= OrderTakeProfit())
    {
         sell_ticket = -1;
    } 
    
    if (OrderSelect(buy_ticket,SELECT_BY_TICKET)==true && OrderStopLoss() != 0 && OrderClosePrice() <= OrderStopLoss())
    {
         buy_ticket = -1;
    } 
    
    if (OrderSelect(sell_ticket,SELECT_BY_TICKET)==true && OrderStopLoss() != 0 && OrderClosePrice() >= OrderStopLoss())
    {
         sell_ticket = -1;
    } 
    
    */
    double ANN_Value = ANN();
   
    if (ANN_Value==1 && buy_ticket == -1)
    {  
           
          
               double TP = 0;
               
               if (TakeProfit != 0)
               {
                  TP = Ask+(TakeProfit*Point*10);
               }
               double SL =  0;
               
               if (StopLoss != 0)
               {
                   SL = Ask-(StopLoss*Point*10);
               }
               
               buy_ticket = OrderSend(Symbol(), OP_BUY, Lot, Ask, OpenSlippage, SL, TP, "Buy Order", Magic, 0, Green);
               
               if (buy_ticket<0)
                {
                A = GetLastError();
                Alert("Buy order opening failed with error: ",ErrorDescription(A));
                A = 0;
                buy_ticket = -1;
                }   
           
             
    }
    
    if (ANN_Value==2 && sell_ticket == -1)
    {
            
           
               TP = 0;
               
               if (TakeProfit != 0)
               {
                  TP = Bid-(TakeProfit*Point*10);
               }
               
               SL =  0;
               
               if (StopLoss != 0)
               {
                   SL = Bid+(StopLoss*Point*10);
               }
         
               sell_ticket = OrderSend(Symbol(), OP_SELL, Lot, Bid, OpenSlippage, SL, TP, "Sell Order", Magic, 0, Red);
               if (sell_ticket<0)
               {
                  A = GetLastError();
                  Alert("Sell order opening failed with error: ",ErrorDescription(A));
                  A = 0;
                  sell_ticket = -1;

               }
        
    
    }
    

    return(0);
 }