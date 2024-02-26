using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;


namespace App
{
    class Description_of_point
    {
        public double x;
        public double y;
        public double Vx;
        public double Vy;

        public Description_of_point(double x, double y, double Vx, double Vy)
        {
            this.x = x;
            this.y = y;
            this.Vx = Vx;
            this.Vy = Vy;
        }
    }
    class Parabola_Flight
    {
        readonly double V0;
        readonly double alpha;
        readonly double m;
        readonly double g = 9.8;
        public List<Description_of_point> Points = new List<Description_of_point> { };
        double partition = 1;
        double k = 1 / 2;

        public Parabola_Flight(double V0, double alpha, double m, double partition)
        {
            this.V0 = V0;
            this.alpha = alpha;
            this.m = m;
            this.partition = partition;
        }

        public Description_of_point Position(List<Description_of_point> Points, Int32 i, double dt,  double m, double k)
        {
            double x = Points[i - 1].x + dt * Points[i - 1].Vx;
            double y = Points[i - 1].y + dt * Points[i - 1].Vy;

            double Vx = Points[i - 1].Vx - k * dt * Points[i - 1].Vx / m;
            double Vy = Points[i - 1].Vy - dt * (g + k * Points[i - 1].Vy / m);

            return new Description_of_point(x, y, Vx, Vy);
        }

        public void Calculate(double x0, double y0)
        {
            double t = 0;
            double T = 2 * V0 * Math.Sin(alpha) / g;
            double dt = T / partition;

            Description_of_point point = new Description_of_point(x0, y0, V0 * Math.Cos(alpha), V0 * Math.Sin(alpha));
            Points.Add(point);
            t += dt;
            Int32 i = 1;

            while (t < T)
            {
                double k = 1 / 2;
                Description_of_point position = Position(Points, i, dt, m, k);
                Points.Add(position);
                t += dt;
                i += 1;
            }
        }
    }

}
