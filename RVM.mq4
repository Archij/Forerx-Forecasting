
//+------------------------------------------------------------------+
//|                                     Relevance Vector Machine.mq4 |
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
extern int Standartization = 0; //0 - ZScore
extern int Normalization = 0; //Нормализация: 0 - [0;1], 1 - [-1;1]
double Average;
double Std;
double x_min;
double x_max;
double candle = 1;

//+------------------------------------------------------------------+
//| expert initialization function                                   |
//+------------------------------------------------------------------+
int buy_ticket = -1;
int sell_ticket = -1;
int error = 0;
bool order_Close = false;
static int prevtime = 0;
int index = 1;
double close[][1];
double Log_Returns[][1];
double Price[][1];
double WMP[][1];
int Training_Data = 261;
bool Train = true;
int embedded_dimension;
int Samples;
double y;
double y_prev;
double DS = 0;
int ds;
double R;
double Covariance = 0;
double STD_Observed;
double STD_Predicted = 0;
double Mean_Observed;
double Mean_Predicted;
double R2;
double SSres = 0;
double SStot = 0;
double RMSE = 0;
double MSE = 0;
double ME = 0;
double RME = 0;
double MAE = 0;
double RMAE = 0;
double MAPE = 0;
 
 
void Transform(double &Time_Series[][], double &Source[][])
{
     int N = ArrayRange(Source,0);
     int indx = 0;
     
     if (Transformation == 0)
     {
        for (int i = 1; i<N; i++)
        {
             Time_Series[indx][0] = MathLog(Source[i][0]/Source[i-1][0]);
             indx++;
        }
     }
     else if (Transformation == 1)
     {
        for (i = 1; i<N; i++)
        {
             Time_Series[indx][0] = Source[i][0]-Source[i-1][0];
             indx++;
        }
     }
     else if (Transformation == 2)
     {
        for (i = 1; i<N; i++)
        {
             Time_Series[indx][0] = Source[i][0]/Source[i-1][0];
             indx++;
        }
     }
}

double Mean(double &Time_Series[][])
{
    int N = ArrayRange(Time_Series,0);
     
    double M = 0;
    
    for (int i=0;i<N;i++)
    {
         M = M + Time_Series[i][0];

    }
    
    M = M/N;
    
    return M;
        
}

double Standart_Deviation(double &Time_Series[][])
{
    int N = ArrayRange(Time_Series,0);
     
    double std = 0;
    
    double mean = Mean(Time_Series);
    
    for (int i=0;i<N;i++)
    {
         std = std + MathPow(Time_Series[i][0]-mean,2);

    }
    
    std = std/(N-1);
    
    std = MathSqrt(std);
    
    return std;
        
}
 
void ZScore(double &Time_Series[][])
{
    
    int N = ArrayRange(Time_Series,0);
    
    double mean = Mean(Time_Series);
    
    double std = Standart_Deviation(Time_Series);
    
    for (int i=0;i<N;i++)
    {
        Time_Series[i][0] = (Time_Series[i][0] - mean) / std;

    }
     
}

void Normalize(double &Time_Series[][], double xmin, double xmax)
{

        int N = ArrayRange(Time_Series,0);
        double a,b;
        
        if (Normalization == 0)
        {
        
            a = 0;
            b = 1;
        
        }
        else if (Normalization == 1)
        {
           a = -1;
           b = 1;    
        }
        
        for (int i=0;i<N;i++)
        {
                Time_Series[i][0]=(((Time_Series[i][0]-xmin)*(b-a))/(xmax-xmin))+a;

        }

}

double Log2(double number)
{
    double log2 = MathLog10(number) / MathLog10(2);
    return(log2);
}

int AMI(double &Time_Series[][], int &tau_values[])
{
   int N = ArrayRange(Time_Series, 0);
   
   double PAB;
   double PA = 0;
   double PB = 0;
   
   double I[][1];
   ArrayResize(I,ArraySize(tau_values),0);
   ArrayInitialize(I,0.0);
   
   for (int i = 0; i<ArraySize(tau_values); i++)
   {
       double x[];
       ArrayResize(x,N-tau_values[i],0);
       
       for (int xi = 0; xi<N-tau_values[i]; xi++)
       {
           x[xi] = Time_Series[xi][0];
       }
           
       double tau[];
       ArrayResize(tau,N-tau_values[i],0);
       
        for (int ti = 0; ti<N-tau_values[i]; ti++)
       {
           tau[ti] = Time_Series[ti+tau_values[i]][0];
       }

       for (int t = 0; t<N-tau_values[i]; t++)
       {
           
           for (int p = 0; p<N-tau_values[i]; p++)
           {
               if (x[t] == x[p])
               {
                   PA = PA + 1;
               }
               
               if (tau[t] == tau[p])
               {
                   PB = PB + 1;
               }
           }
           

           PAB = (PA + PB)/ (2*(N-tau_values[i]));
           PA = PA/(N-tau_values[i]);
           PB = PB/(N-tau_values[i]);

           I[i][0] = I[i][0] + PAB * Log2(PAB/(PA*PB));
           
           PA = 0;
           PB = 0;
      
          
       }
   }
   
   int min = ArrayRange(I,0)-1;
   
   for (int m = 1; m<ArraySize(tau_values); m++)
   {
       if (I[m][0] > I[m-1][0])
       {
          min = m-1;
          break;
       }
   }
   
   return(min);
}

int FNN(double &Time_Series[][], int tau)
{

    
    
}

void TD_PhaseSpace(double &Time_Series[][],int m,int tau,double &Y[][])
{
    int N=ArrayRange(Time_Series,0);

    int M=N-(m-1)*tau;

    for (int i=0; i<m; i++)
    {
        for (int j=0; j<M; j++)
        {
            Y[j][i] = Time_Series[j+i*tau][0];               
        }
    }
}


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
  int n = ArrayRange(Matrix1, 1);
  int m = ArrayRange(Matrix1, 0);
  int p = ArrayRange(Matrix2, 1);
  
     for (int i=0; i<m; i++)
     {
         for (int j=0; j<p; j++)
         {
             for (int k=0; k<n; k++)
             {
                 Result[i][j] = Result[i][j] + Matrix1[i][k] * Matrix2[k][j];
             }
         }
     }

}

void ScalarMultiply(double value, double &Matrix[][], double &Result[][])
{
    for (int i = 0; i<ArrayRange(Matrix, 0); i++)
    {
        for (int j = 0; j<ArrayRange(Matrix, 1); j++)
        {
            
               Result[i][j] = Matrix[i][j]*value;
            
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

void add(double &Matrix1[][], double &Matrix2[][], double &Result[][])
{

       for (int i = 0; i<ArrayRange(Matrix1, 0); i++)
       {
           for (int j = 0; j<ArrayRange(Matrix1, 1); j++)
           {
               Result[i][j] = Matrix1[i][j] + Matrix2[i][j];
           }
       }   
       
}

void subtract(double &Matrix1[][], double &Matrix2[][], double &Result[][])
{
       for (int i = 0; i<ArrayRange(Matrix1, 0); i++)
       {
           for (int j = 0; j<ArrayRange(Matrix1, 1); j++)
           {
               Result[i][j] = Matrix1[i][j] - Matrix2[i][j];
           }
       }
        
}

void diag(double &Vector[][], double &Matrix[][])
{
    for (int i = 0; i<ArrayRange(Matrix, 0); i++)
    {
        for (int j = 0; j<ArrayRange(Matrix, 1); j++)
        {
            if (i == j)
            {
               Matrix[i][j] = Vector[i][0];
            }
        }
    }
}

double Norm(double &Vector[][])
{
    double norm = 0;
    
    for (int v = 0; v < ArrayRange(Vector,0); v++)
    {
       norm = norm + MathPow(MathAbs(Vector[v][0]),2);
    }
    
    return(MathSqrt(norm));
}

double K(double &X[][],double &Y[][],string kernel_)
{

    double f = 0;
    int N = ArrayRange(X,0);
    double c = 0;
    double XT[][2];
    ArrayResize(XT,1,0);
    double Result[][1];
    ArrayResize(Result,1,0);
    double alpha = 1;
    int d = 2;
    double q = 0.1;
    double matrix[][1];
    ArrayResize(matrix,N,0);
    double norm =  0;
    int n = N;
    double Q = 1;
    double C = 1;
    double a = 1;
    
    if (kernel_ == "Linear")
    {
           ArrayInitialize(XT,0.0);
           Transpose(X,XT);
           ArrayInitialize(Result,0.0);
           Multiply(XT,Y,Result);
           f = Result[0][0]+c;
    }
    else if (kernel_ == "Polynomial")
    {
           ArrayInitialize(XT,0.0);
           Transpose(X,XT);
           ArrayInitialize(Result,0.0);
           Multiply(XT,Y,Result);
           f = MathPow(alpha*Result[0][0]+c,d);
    }
    else if (kernel_ == "Gauss")
    {
           subtract(X,Y,matrix);
           norm = MathPow(Norm(matrix),2);
           f = MathExp(-(norm / (2*MathPow(q,2))));
    }
    else if (kernel_ == "Exponential")
    {
           subtract(X,Y,matrix);
           norm = Norm(matrix);
           f = MathExp(-(norm / (2*MathPow(q,2))));
    }
    else if (kernel_ == "Laplacian")
    {
           subtract(X,Y,matrix);
           norm = Norm(matrix);
           f = MathExp(-(norm / q));
    }
    else if (kernel_ == "Sigmoid")
    {
           ArrayInitialize(XT,0.0);
           Transpose(X,XT);
           ArrayInitialize(Result,0.0);
           Multiply(XT,Y,Result);
           f = MathTan(alpha*Result[0][0]+c);
    }
    else if (kernel_ == "Rational_Quadratic")
    {
           subtract(X,Y,matrix);
           norm = MathPow(Norm(matrix),2);
           f = 1 - norm/(norm + c);
    }
    else if (kernel_ == "Multiquadric")
    {
           subtract(X,Y,matrix);
           norm = MathPow(Norm(matrix),2);
           f = MathSqrt(norm + MathPow(c,2));
    }
    else if (kernel_ == "Inverse_Multiquadric")
    {
           subtract(X,Y,matrix);
           norm = MathPow(Norm(matrix),2);
           f = 1/MathSqrt(norm + MathPow(c,2));
    }
    else if (kernel_ == "Circular")
    {
           subtract(X,Y,matrix);
           norm = Norm(matrix);
           if (norm < q)
           {
              f = (2/M_PI)*MathArccos(-norm/q)-(2/M_PI)*(norm/q)*MathSqrt(1 - MathPow(norm/q,2));
           }  
    }
    else if (kernel_ == "Spherical")
    {
           subtract(X,Y,matrix);
           norm = Norm(matrix);
           if (norm < q)
           {
              f = 1 - (3/2)*(norm/q) + (1/2)*MathPow(norm/q,3);
           }
    }
    else if (kernel_ == "Wave")
    {
           subtract(X,Y,matrix);
           norm = Norm(matrix);
           f = (Q/norm)*MathSin(norm/Q);  
    }
    else if (kernel_ == "Power")
    {
           subtract(X,Y,matrix);
           norm = Norm(matrix);
           f = -MathPow(norm,d);     
    }
    else if (kernel_ == "Log")
    {
           subtract(X,Y,matrix);
           norm = Norm(matrix);
           f = -MathLog(MathPow(norm,d)+1);     
    }
    else if (kernel_ == "Cauchy")
    {
           subtract(X,Y,matrix);
           norm = MathPow(Norm(matrix),2);
           f = 1 / (1 + norm/MathPow(q,2));     
    }
    else if (kernel_ == "Chi-Square")
    {
           for (int i = 0; i<n; i++)
           {
                f = f + MathPow(X[i][0]-Y[i][0],2)/((1/2)*(X[i][0]+Y[i][0]));
           }
           
           f = 1 - f;     
    }
    else if (kernel_ == "T-Student")
    {
           subtract(X,Y,matrix);
           norm = MathPow(Norm(matrix),d);
           f = 1 / (1 + norm);     
    }
    else if (kernel_ == "Wavelet")
    {
           f = 1;
           for (i = 0; i<N; i++)
           {
               f = f * ((MathCos(1.75*((X[i][0]-C)/a))*MathExp(-MathPow((X[i][0]-C)/a,2)/2)) * (MathCos(1.75*((Y[i][0]-C)/a))*MathExp(-MathPow((Y[i][0]-C)/a,2)/2)));
           }   
    }
    else if (kernel_ == "Translation-Invariant_Wavelet")
    {
           f = 1;
           for (i = 0; i<N; i++)
           {
               f = f * (MathCos(1.75*((X[i][0]-Y[i][0])/a))*MathExp(-MathPow((X[i][0]-Y[i][0])/a,2)/2));
           }   
    }
    
    return(f);
}

void PhiMatrix(double &X[][], double &Phi[][])
{
     int N = ArrayRange(Phi,0);
     int M = ArrayRange(Phi,1);
     int d = ArrayRange(X,1);
     
     double X1[][1];
     ArrayResize(X1,d,0);
     double X2[][1];
     ArrayResize(X2,d,0);
     
     
     for (int n = 0; n<N; n++)
     {
         for (int m = 0; m<M; m++)
         {
             if (m == 0)
             {
                Phi[n][m] = 1; //bias
             }
             else
             {
                for (int ki = 0; ki<d; ki++)
                {
                    X1[ki][0] = X[n][ki];
                    X2[ki][0] = X[m-1][ki];
                }
                
                Phi[n][m] = K(X1,X2,"Gauss");
             }
         }
     }
}

void RVR(double &Phi[][], double &t[][])
{
      // Terminate estimation when no log-alpha value changes by more than this
	   const double Delta_Min = 1e-6;

	   // Prune basis function when its alpha is greater than this
	   const double Alpha_Max = 1e12;
	   
      int N = ArrayRange(Phi, 0);
      int M = ArrayRange(Phi, 1);
	  
	   int maxIts = 1000;
	   
	   double initAlpha = MathPow(1.0 / double(N),2);
	   double epsilon = Standart_Deviation(t)*10.0/100.0;
	   double Beta = 1.0 / MathPow(epsilon,2);
      
      double Alpha[][1];
      ArrayResize(Alpha,M);

      ArrayResize(WMP,M,0);
      
      for (int w = 0; w<M; w++)
      {
         WMP[w][0] = 0;
      }
      
      for (int a = 0; a<M; a++)
      {
         Alpha[a][0] = initAlpha;
      }
      
      double Gamma[][1];
      ArrayResize(Gamma, M, 0);
      
      for (int g = 0; g<M; g++)
      {
         Gamma[g][0] = 1.0;
      }
      
      double Delta[][1];
      ArrayResize(Delta, M, 0);

      double PhiT[][131]; 
      ArrayResize(PhiT,M,0);
      
      Transpose(Phi, PhiT);
      
      for (int k = 1; k<=maxIts; k++)
      {
          Alert("iteration: ", k);
          double A[][132];
          ArrayResize(A,M,0);
          ArrayInitialize(A,0.0);
          
          diag(Alpha, A);
          
           //Alert("A ", A[0][0], " ", A[0][1], " ", A[0][2], " ", A[1][0], " ", A[1][1], " ", A[1][2], " ", A[2][0], " ", A[2][1], " ", A[2][2]);

          double Sigma[][132];
          ArrayResize(Sigma,M,0);
          
          double bPhiT[][131];
          ArrayResize(bPhiT,M,0);
          
          ScalarMultiply(Beta, PhiT, bPhiT);
          
          double bPhiTPhi[][132];
          ArrayResize(bPhiTPhi,M,0);
          ArrayInitialize(bPhiTPhi,0.0);
          
          Multiply(bPhiT, Phi, bPhiTPhi);
          
          double bPhiTPhiA[][132];
          ArrayResize(bPhiTPhiA,M,0);
          
          add(bPhiTPhi, A, bPhiTPhiA);

          Inverse(bPhiTPhiA, Sigma);
          
          //Alert("Sigma ", Sigma[0][0], " ", Sigma[0][1], " ", Sigma[0][2], " ", Sigma[1][0], " ", Sigma[1][1], " ", Sigma[1][2], " ", Sigma[2][0], " ", Sigma[2][1], " ", Sigma[2][2]);
          
          double Sigmab[][132];
          ArrayResize(Sigmab,M,0);
          
          ScalarMultiply(Beta, Sigma, Sigmab);
         
          double SigmabPhiT[][131]; 
          ArrayResize(SigmabPhiT,M,0);
          ArrayInitialize(SigmabPhiT,0.0);
      
          Multiply(Sigmab, PhiT, SigmabPhiT);
          
          ArrayInitialize(WMP,0.0);
          Multiply(SigmabPhiT, t, WMP);
          
          //Alert("WMP", WMP[0][0], " ", WMP[1][0], " ", WMP[2][0]);
          
          double AlphaOld[][1];
          ArrayResize(AlphaOld,M);
          
          for (int d = 0; d<M; d++)
          {
              Delta[d][0] = 0;
          }
          
          for (int j = 0; j<M; j++)
          {
       
              if (Alpha[j][0] > Alpha_Max)
              {
                  WMP[j][0] = 0;
                  Alpha[j][0] = MathExp(1000);
                  Gamma[j][0] = 0;
              }
              else
              {   
                  AlphaOld[j][0] = Alpha[j][0];
                  Gamma[j][0] = 1 - Alpha[j][0] * Sigma[j][j];
                  if (k < maxIts/2)
                  {
                       Alpha[j][0] = Gamma[j][0] / MathPow(WMP[j][0],2);
                  }
                  else
                  {
                       Alpha[j][0] = Gamma[j][0] * MathPow(MathPow(WMP[j][0],2)/Gamma[j][0]-Sigma[j][j],-1);
                  }
                  Delta[j][0] = MathAbs(MathLog(Alpha[j][0]) - MathLog(AlphaOld[j][0]));
              }
              
          }
          
          //Alert("WMP ", WMP[0][0], " ", WMP[1][0], " ", WMP[2][0]);
          //Alert("Alpa ", Alpha[0][0], " ", Alpha[1][0], " ", Alpha[2][0]);
          //Alert("Gamma ", Gamma[0][0], " ", Gamma[1][0], " ", Gamma[2][0]);
          
          
          double gsum = 0;
          for (j = 0; j<M; j++)
          {
              gsum = gsum + Gamma[j][0];
          }
          
          double PhiWMP[][1];
          ArrayResize(PhiWMP,N,0);
          ArrayInitialize(PhiWMP,0.0);
          
          Multiply(Phi, WMP, PhiWMP);
          
          double tPhiWMP[][1];
          ArrayResize(tPhiWMP,N,0);
          
          subtract(t, PhiWMP, tPhiWMP);
          
          double norm = Norm(tPhiWMP); 

          Beta = (N - gsum) / MathPow(norm,2);

          if (Delta[ArrayMaximum(Delta,WHOLE_ARRAY)][0] < Delta_Min)
          {
             break;

          }

      }
      

       //Alert("iteration: ", k-1);
       //Alert("WMP: ", WMP[0][0], " ", WMP[1][0], " ", WMP[2][0], " ", WMP[3][0], " ", WMP[4][0], " ", WMP[5][0]);
       //Alert("Y: ", PhiWMP[0][0], " ", PhiWMP[1][0], " ", PhiWMP[2][0], " ", PhiWMP[3][0], " ", PhiWMP[4][0], " ");
      // Alert("t: ", t[0][0], " ", t[1][0], " ", t[2][0], " ", t[3][0], " ", t[4][0], " ");
       

       Mean_Observed = Mean(t);
       Mean_Predicted = Mean(PhiWMP);
       
       for (int ti=0; ti<N; ti++)
       {
           if (ti > 0)
           {
              if (((PhiWMP[ti][0] - PhiWMP[ti-1][0]) * (t[ti][0] - t[ti-1][0]))>0)
              {
                 ds = 1;
              }
              else
              {
                 ds = 0;
              }
           }
           
           DS = DS + ds;
           Covariance = Covariance + (t[ti][0]-Mean_Observed)*(PhiWMP[ti][0]-Mean_Predicted);
           STD_Observed = STD_Observed + MathPow(t[ti][0]-Mean_Observed,2);
           STD_Predicted = STD_Predicted + MathPow(PhiWMP[ti][0]-Mean_Predicted,2);
           SSres = SSres + MathPow(t[ti][0]-PhiWMP[ti][0],2);
           SStot = SStot + MathPow(t[ti][0]-Mean_Observed,2);
           MSE = MSE + MathPow(t[ti][0]-PhiWMP[ti][0],2);
           ME = ME + (t[ti][0]-PhiWMP[ti][0]);
           if (t[ti][0] != 0)
           {
              RME = RME + (PhiWMP[ti][0]-t[ti][0])/t[ti][0];
              RMAE = RMAE + MathAbs(PhiWMP[ti][0]-t[ti][0])/t[ti][0];
              MAPE = MAPE + MathAbs((t[ti][0]-PhiWMP[ti][0])/t[ti][0])*100/N;
           }
           else
           {
              RME = RME + MathExp(1000);
              RMAE = RMAE + MathExp(1000);
              MAPE = MAPE + MathExp(1000);
           }
           MAE = MAE + MathAbs(t[ti][0]-PhiWMP[ti][0]);
           

           
       }
       
       DS = DS*100/(N-1);
       if (MathSqrt(STD_Observed)*MathSqrt(STD_Predicted) != 0)
       {
          R = (Covariance) / (MathSqrt(STD_Observed)*MathSqrt(STD_Predicted));
       }
       else
       {
          R = MathExp(1000);
       }
       if (SStot != 0)
       {
          R2 = 1 - SSres/SStot;
       }
       else
       {
          R2 = MathExp(1000);
       }
       MSE = MSE/N;
       RMSE = MathSqrt(MSE);
       RME = RME/N;
       MAE = MAE/N;
       
       Alert("DS: ", DS);
       Alert("R: ", R);
       Alert("R^2: ", R2);
       Alert("RMSE: ", RMSE);
       Alert("RME: ", RME);
       Alert("MAE: ", MAE);
       Alert("MAPE: ", MAPE);
       
       DS = 0;
            
}

int init()
  {

     prevtime = iTime(Symbol(),0,0);
     
    
     ArrayResize(close,Training_Data,0);
    
    for (int i = 0; i<Training_Data; i++)
    {
        close[i][0] = Close[Training_Data-i];
    }
    
    Alert(DoubleToStr(close[0][0],5));
    Alert(DoubleToStr(close[260][0],5));

    ArrayResize(Price, Training_Data-1, 0);

    //Трансформация
    Transform(Price, close);
    
    //Стандартизация
    //Average = Mean(Price);
    //Std = Standart_Deviation(Price);
    //ZScore(Price);
    
    //Нормализация
    //x_min=Price[ArrayMinimum(Price,WHOLE_ARRAY,0)][0];
    //x_max=Price[ArrayMaximum(Price,WHOLE_ARRAY,0)][0];
    //Normalize(Price, x_min, x_max);
    
    int N = ArrayRange(Price,0);
    
    embedded_dimension = 129;
    
    Samples = N - embedded_dimension;
    
    double X[][129];
    ArrayResize(X,Samples,0);
       
    double t[][1];
    ArrayResize(t,Samples,0);
    
    
    for (int s = 0; s<Samples; s++)
    {
        for (int d = 0; d<embedded_dimension; d++)
        {
            X[s][d] = Price[s+d][0];
        }
        
        t[s][0] = Price[s+embedded_dimension][0];
    }
         
    double PHI[][132];
    ArrayResize(PHI,Samples,0);
    
    PhiMatrix(X, PHI);
    
    int M = ArrayRange(PHI,1);
    ArrayResize(WMP, M, 0);
  
    RVR(PHI,t);

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

    for (int i = 0; i<Training_Data; i++)
    {
        close[i][0] = Close[Training_Data-i];
    }

    //Трансформация
    Transform(Price, close);
    
    //Average = Mean(Price);
    //Std = Standart_Deviation(Price);
    //ZScore(Price);
    
        //Нормализация
    //x_min=Price[ArrayMinimum(Price,WHOLE_ARRAY,0)][0];
    //x_max=Price[ArrayMaximum(Price,WHOLE_ARRAY,0)][0];
    //Normalize(Price, x_min, x_max);
    
    
    double X[][129];
    ArrayResize(X,Samples,0);
    
    for (int s = 0; s<Samples; s++)
    {
        for (int d = 0; d<embedded_dimension; d++)
        {
            X[s][d] = Price[s+d+1][0];
        }
    }
    
    /*
    double K[][1];
    ArrayResize(K,Samples,0);
    
    for (int k = 0; k<Samples; k++)
    {
        for (int ki = 0; ki<embedded_dimension; ki++)
        {
             X1[ki][0] = X[n][ki];
             X2[ki][0] = X[m-1][ki];
        }
                
         K[k][0] = K(X1,X2,"Gauss");
    }
    */
    double PHI[][132];
    ArrayResize(PHI,Samples,0);
         
    PhiMatrix(X, PHI);

    double PhiWMP[][1];
    ArrayResize(PhiWMP,Samples,0);
    ArrayInitialize(PhiWMP,0.0);
          
    Multiply(PHI, WMP, PhiWMP);
    
    /*
    double a = x_min;
    double b = x_max;
        
    if (Normalization == 0)
    {
            x_min = 0;
            x_max = 1;
    }
    else if (Normalization == 1)
    {
           x_min = -1;
           x_max = 1;    
    }

    y = (((PhiWMP[Samples-1][0]-x_min)*(b-a))/(x_max-x_min))+a;
    
    y = y*Std + Average;
    
   
    */
    
    y = MathExp(PhiWMP[Samples-1][0])*Close[1];
    
    if (buy_ticket != -1)
    {
            
              order_Close = OrderClose(buy_ticket,Lot,Bid,CloseSlippage,Green);
           
              if (order_Close == false)
              {
                  error = GetLastError();
                  Alert("Buy order closing failed with error: ",ErrorDescription(error));
                  error = 0;
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
                    error = GetLastError();
                    Alert("Sell order closing failed with error: ",ErrorDescription(error));
                    error = 0;
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
    


   
    if (y > Open[0] && buy_ticket == -1)
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
               /*
               Lot = (Risk*AccountFreeMargin())/100 * AccountLeverage() * Close[0] / 100000;
            
           
                  if (Lot > MarketInfo(Symbol(), MODE_MAXLOT))
                  {
                     Lot = MarketInfo(Symbol(), MODE_MAXLOT);
                  }
                  else if (Lot < MarketInfo(Symbol(), MODE_MINLOT))
                  {
                     Lot = MarketInfo(Symbol(), MODE_MINLOT);
                  }
               */
               buy_ticket = OrderSend(Symbol(), OP_BUY, Lot, Ask, OpenSlippage, SL, TP, "Buy Order", Magic, 0, Green);
               
               if (buy_ticket<0)
                {
                error = GetLastError();
                Alert("Buy order opening failed with error: ",ErrorDescription(error));
                error = 0;
                buy_ticket = -1;
                }   
           
             
    }
    
    if (y < Open[0] && sell_ticket == -1)
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
               /*
                Lot = (Risk*AccountFreeMargin())/100 * AccountLeverage() * Close[0] / 100000;
            
           
                  if (Lot > MarketInfo(Symbol(), MODE_MAXLOT))
                  {
                     Lot = MarketInfo(Symbol(), MODE_MAXLOT);
                  }
                  else if (Lot < MarketInfo(Symbol(), MODE_MINLOT))
                  {
                     Lot = MarketInfo(Symbol(), MODE_MINLOT);
                  }
         */
               sell_ticket = OrderSend(Symbol(), OP_SELL, Lot, Bid, OpenSlippage, SL, TP, "Sell Order", Magic, 0, Red);
               if (sell_ticket<0)
               {
                  error = GetLastError();
                  Alert("Sell order opening failed with error: ",ErrorDescription(error));
                  error = 0;
                  sell_ticket = -1;

               }
        
    
    }
    
    candle = candle + 1;
    
    return(0);
 }