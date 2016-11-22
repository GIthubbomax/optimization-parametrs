using System;
using System.IO;
using System.Threading;

namespace optimization_parametrs
{
    class Program
    {
        static void Main(string[] args)
        {
            var a = new {nan= "нан", jopa = true};
            Console.WriteLine(a.jopa+a.nan);
            File.AppendAllText("C:\\Users\\Бочков\\Desktop\\results\\restab.csv",
                "thread" + ";" + "keylength" + ";" + "mu" + ";" + "nu" + ";" + "pmu" + ";" + "pnu" + ";" +
                "plambda" + ";" + "input_mu" + ";" + "input_nu" + ";" + "input_pmu" + ";" + "input_pnu" + ";" +
                "input_plambda" + System.Environment.NewLine);
            for (int i = 0; i < 2; i++)
            {
                KeyLength keyLength = new KeyLength(i.ToString());
                Thread th = new Thread(keyLength.Run);
                th.Start();
                Thread.Sleep(100);
            }


            Console.ReadLine();
        }
    }

    class Answer
    {
        public Answer(double f, double mu, double nu, double pmu, double pnu, double plambda)
        {
            this.mu = mu;
            this.f = f;
            this.nu = nu;
            this.pmu = pmu;
            this.pnu = pnu;
            this.plambda = plambda;
        }


        public double mu { get; set; }
        public double nu { get; set; }
        public double pmu { get; set; }
        public double pnu { get; set; }
        public double plambda { get; set; }

        public double f { get; set; }

        public override string ToString()
        {
            return "Длина ключа: " + String.Format("{0:0.00}", f) + " mu " + String.Format("{0:0.00}", mu) + " nu " +
                   String.Format("{0:0.00}", nu) + " p_mu " + String.Format("{0:0.00}", pmu) + " p_nu " +
                   String.Format("{0:0.00}", pnu) + " p_lambda " +
                   String.Format("{0:0.00}", plambda);
        }
    }

    class KeyLength
    {
        object lockOn = new object();
        private Answer answer;


        private string thteadname;
        double mu, nu, pnu, pmu, plambda;

        public KeyLength(string s)
        {
            thteadname = s;
        }

        const double eta = 0.0055603;
        const double pdc = 0.000062;
        const double popt = 0.025;

        const double mumax = 0.9999999999999;

        const double eps = 0.01;
        const double R = 0.9;
        const double EPSpa = 0.000000000001;
        const double lambda = 0.01;
        const double N = 1372669434;
        const double fi = 14.764;
        const double n = 995328;

        double N_nu(double p_nu)
        {
            return p_nu*N;
        }

        double N_mu(double p_mu)
        {
            return p_mu*N;
        }

        double N_lambda(double p_lambda)
        {
            return p_lambda*N;
        }

        double Q_nu(double nu)
        {
            return nu*eta + 2*pdc;
        }

        double Q_lambda()
        {
            return lambda*eta + 2*pdc;
        }

        double Q_mu(double mu)
        {
            return mu*eta + 2*pdc;
        }

        double Q_Unu(double nu, double p_nu)
        {
            return Q_nu(nu) + fi*Math.Sqrt((Q_nu(nu)*(1 - Q_nu(nu)))/N_nu(p_nu));
        }

        double Q_Umu(double mu, double p_mu)
        {
            return Q_mu(mu) + fi*Math.Sqrt((Q_mu(mu)*(1 - Q_mu(mu)))/N_mu(p_mu));
        }

        double Q_Ulambda(double p_lambda)
        {
            return Q_lambda() + fi*Math.Sqrt((Q_lambda()*(1 - Q_lambda()))/N_lambda(p_lambda));
        }

        double Q_Lnu(double nu, double p_nu)
        {
            return Q_nu(nu) - fi*Math.Sqrt((Q_nu(nu)*(1 - Q_nu(nu)))/N_nu(p_nu));
        }

        double Q_Lmu(double mu, double p_mu)
        {
            return Q_mu(mu) - fi*Math.Sqrt((Q_mu(mu)*(1 - Q_mu(mu)))/N_mu(p_mu));
        }

        double Q_Llambda(double p_lambda)
        {
            return Q_lambda() - fi*Math.Sqrt((Q_lambda()*(1 - Q_lambda()))/N_lambda(p_lambda));
        }

        double Y_Lo(double nu, double p_lambda, double p_nu)
        {
            return Math.Max(0,
                (nu*Q_Llambda(p_lambda)*Math.Exp(lambda) - lambda*Q_Unu(nu, p_nu)*Math.Exp(nu))/(nu - lambda));
        }

        private double hx(double x)
        {
            double a = -x*Math.Log(x, 2) - (1 - x)*Math.Log(1 - x, 2);
            //   Console.WriteLine("hx {0}", a);
            return a;
        }


        private double e_mu(double mu)
        {
            double a = (mu*eta*popt + pdc)/(mu*eta + 2*pdc);
            // Console.WriteLine("функция e_mu :{0}", a);
            return a;
        }

        double Q_L1(double mu, double nu, double p_nu, double p_lambda, double p_mu)
        {
            return mu*Math.Exp(-mu)/(nu*(1 - nu/mu) - lambda*(1 - lambda/mu))*
                   (Q_Lnu(nu, p_nu)*Math.Exp(nu) - Q_Ulambda(p_lambda)*Math.Exp(lambda) -
                    (nu*nu - lambda*lambda)/(mu*mu)*(Q_Umu(mu, p_mu)*Math.Exp(mu) - Y_Lo(nu, p_lambda, p_nu)));
        }

        double E_Umu(double mu)
        {
            return e_mu(mu) + fi*Math.Sqrt(e_mu(mu)*(1 - e_mu(mu))/n);
        }

        double tetha_L1(double mu, double nu, double p_nu, double p_lambda, double p_mu)
        {
            return Q_L1(mu, nu, p_nu, p_lambda, p_mu)/Q_Umu(mu, p_mu);
        }

        double E_U1(double mu, double nu, double p_nu, double p_lambda, double p_mu)
        {
            return (E_Umu(mu)*Q_Umu(mu, p_mu) - Y_Lo(nu, p_lambda, p_nu)*Math.Exp(-mu)/2)/
                   Q_L1(mu, nu, p_nu, p_lambda, p_mu);
        }

        double K_L1(double mu, double nu, double p_nu, double p_lambda, double p_mu)
        {
            return tetha_L1(mu, nu, p_nu, p_lambda, p_mu) -
                   fi*Math.Sqrt(tetha_L1(mu, nu, p_nu, p_lambda, p_mu)*(1 - tetha_L1(mu, nu, p_nu, p_lambda, p_mu))/n);
        }

        double e_U1(double mu, double nu, double p_nu, double p_lambda, double p_mu)
        {
            return E_U1(mu, nu, p_nu, p_lambda, p_mu) +
                   fi*Math.Sqrt(E_U1(mu, nu, p_nu, p_lambda, p_mu)*(1 - E_U1(mu, nu, p_nu, p_lambda, p_mu))/n);
        }

        public double l(double mu, double nu, double p_mu, double p_nu, double p_lambda)
        {
            return K_L1(mu, nu, p_nu, p_lambda, p_mu)*n*(1 - hx(e_U1(mu, nu, p_nu, p_lambda, p_mu))) - 1.22*hx(e_mu(mu)) -
                   5*Math.Log(1/EPSpa, 2);
        }

        public void Run()
        {
            Random random = new Random();

            mu = random.NextDouble();
            nu = random.NextDouble();
            pnu = pmu = plambda = -1;
            while (pmu <= 0 || pmu > mumax || pnu <= 0 || pnu > mumax || plambda < 0 ||
                   plambda > mumax ||
                   (pnu + pmu + plambda) !=1.000000000000000000)
            {
                pnu = random.NextDouble();
                pmu = random.NextDouble();
                plambda = 1 - pnu - pmu;
            }
            Console.WriteLine("данные поток {5} mu={0:0.00},nu={1:0.00},pmu={2:0.00},pnu={3:0.00},plambda={4:0.00}", mu,
                nu, pmu, pnu, plambda, thteadname);

            answer = annealing(mu, nu, pmu, pmu, plambda);
 //           answer = annealing(0.03,0.012,0.366,0.784,0.555);
            Console.WriteLine("результат поток {0} " + answer, thteadname);
            Write();
        }

        void Write()
        {
            lock (lockOn)
            {
                File.AppendAllText("C:\\Users\\Бочков\\Desktop\\results\\restab.csv",
                    String.Format("{0};{1};{2};{3};{4};{5};{6};{7};{8};{9};{10};{11}" + System.Environment.NewLine,
                        thteadname, answer.f, answer.mu, answer.nu, answer.pmu, answer.pnu, answer.plambda, mu, nu, pmu,
                        pnu, plambda));
            }
        }

        public Answer annealing(double mu, double nu, double pmu, double pnu, double plambda)
        {
            double munew, nunew, pnunew, pmunew, plambdanew, alpha, z;
            int a = 0;
            double h, f, fnew;
            f = l(mu, nu, pmu, pnu, plambda);
            double temperature = 10000;
            while (temperature > eps && a < 10000000000)
            {
                a++;
                Random random = new Random();
                pmunew = plambdanew = pnunew = munew = nunew = -1;

                while (munew <= 0 || munew > mumax)
                {
                    alpha = random.NextDouble();
                    z = Math.Sign(alpha - 0.5)*temperature*(Math.Pow((1 + 1/temperature), Math.Abs(2*alpha - 1)) - 1);
                    munew = mu + z*mumax;
                }
                while (nunew <= 0 || nunew > mumax)
                {
                    alpha = random.NextDouble();
                    z = Math.Sign(alpha - 0.5)*temperature*(Math.Pow((1 + 1/temperature), Math.Abs(2*alpha - 1)) - 1);
                    nunew = nu + z*mumax;
                }
                while (pmunew <= 0 || pmunew > mumax || pnunew <= 0 || pnunew > mumax || plambdanew < 0 ||
                       plambdanew > mumax ||
                       (pnunew + pmunew + plambdanew) != 1.0)
                {
                    alpha = random.NextDouble();
                    z = Math.Sign(alpha - 0.5)*temperature*(Math.Pow((1 + 1/temperature), Math.Abs(2*alpha - 1)) - 1);
                    pnunew = pnu + z*mumax;
                    alpha = random.NextDouble();
                    z = Math.Sign(alpha - 0.5)*temperature*(Math.Pow((1 + 1/temperature), Math.Abs(2*alpha - 1)) - 1);
                    pmunew = pmu + z*mumax;
                    plambdanew = 1 - pmunew - pnunew;
                }


                fnew = l(munew, nunew, pmunew, pnunew, plambdanew);

                if (fnew > f)
                {
                    mu = munew;
                    nu = nunew;
                    pnu = pnunew;
                    pmu = pmunew;
                    plambda = plambdanew;
                    temperature = temperature*R;
                    f = fnew;
                }
                else
                {
                    h = 1/
                        (1 +
                         Math.Exp(-((fnew - f))/
                                  temperature));
                    if (h > random.NextDouble())
                    {
                        mu = munew;
                        nu = nunew;
                        pnu = pnunew;
                        pmu = pmunew;
                        plambda = plambdanew;
                        temperature = temperature*R;
                        f = fnew;
                    }
                }
            }
            return new Answer(l(mu, nu, pmu, pnu, plambda), mu, nu, pmu, pnu, plambda);
        }
    }
}