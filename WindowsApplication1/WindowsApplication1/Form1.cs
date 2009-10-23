using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Collections;
using System.Windows.Forms;

namespace WindowsApplication1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
			getU();
            InitializeComponent();
        }
		/*
        public void SFFT(int Bsample, int l, int m)
        {
            int B2 = 2*Bsample - 1;
            for(int j=0; j < B2; j++)
            {
                double thetaJ = (System.Math.PI*(2*j+1))/(4*Bsample);
                for (int k = 0; k < B2; k++)
                {
                    double phiK = (System.Math.PI * 2 * k) / (2 * B);
                    int a = methode(j,Bsample);
                    Vec3D fSample = new Vec3D(System.Math.Cos(phiK) * System.Math.Sin(thetaJ),
                                                    System.Math.Sin(phiK) * System.Math.Sin(thetaJ),
                                                    System.Math.Cos(thetaJ)     );
                    int Y = aargh();
                }
            }
        }*/

		public void getU()
		{
			
		}
    }
}