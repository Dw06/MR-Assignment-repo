/*
 * Multimedia Retrieval Practical Assignment: Ray-Based Approach with Spherical Harmonics
 * By Milan de Graaf, 3117308 and David Weterings, 3117480
 * 
 **/ 
using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Windows.Forms;
using System.Runtime.InteropServices;

namespace WindowsApplication1
{

    static class Program
    {
        //possible values for B: 32,64,128,256,512,1024 etc.
        const int B = 64;

        const int B2 = B * 2;

        //dll functions
        [DllImport("SphericTrans.dll")]
        public static extern unsafe void calcCoeff(int bw, int cutoff, [MarshalAs(UnmanagedType.LPArray)] double[] rdata, [MarshalAs(UnmanagedType.LPArray)] double[] idata, [MarshalAs(UnmanagedType.LPArray)] double[] rcoeffs, [MarshalAs(UnmanagedType.LPArray)] double[] icoeffs);

        [DllImport("SphericTrans.dll")]
        public static extern unsafe void inverseTransform(int bw, int cutoff, [MarshalAs(UnmanagedType.LPArray)] double[] rcoeffs, [MarshalAs(UnmanagedType.LPArray)] double[] icoeffs, [MarshalAs(UnmanagedType.LPArray)] double[] rdata, [MarshalAs(UnmanagedType.LPArray)] double[] idata);

        //model directory
        const string modelPath = "../../test data/";
		
        /// <summary>
        /// The main entry point for the application.
        /// </summary>
        [STAThread]
        static void Main()
        {
            DateTime startTimeTotal = DateTime.Now;
            Console.WriteLine("parsing basenames");
            ArrayList files = FileParser.parseBaseFileSimple(modelPath, "basenames");
            Console.WriteLine("Get all directions");
            unitVector[] directSorted = getDirections();
            //sort on phi, theta
            Array.Sort(directSorted); 
            Console.WriteLine("Start comparing");

            foreach(String f in files)
            {
                Console.WriteLine("Loading model: " + f);
                int total = 0;
                double[] distances = new double[files.Count];
                Model3D curModel = FileParser.parseModel(modelPath + f, f);
                Console.WriteLine("Extracting fv: " + f);
                if (curModel == null)
                    continue;
                double[] FvCurrent;
                string sourceFileC = modelPath + f + ".fv_B_" + B;
                if (File.Exists(sourceFileC))
                {
                    Console.WriteLine("Obtaining feature vector from file");
                    FvCurrent = FileParser.parseFeatureVector(modelPath, f, B);
                }
                else
                {
                    Console.WriteLine("Calculating feature vector");
                    FvCurrent = extractFeatureVector(ref curModel, ref directSorted);
                    FileWriter.WriteFeatureFile(FvCurrent, modelPath, f);
                }
               
                foreach (String s in files)
                {
                    Console.WriteLine("Loading other model: " + s);
                    Model3D tarModel = FileParser.parseModel(modelPath + s, s);
                    if(tarModel == null)
                    {
                        distances[total] = -1;
                        total++;
                        continue;
                    }
                    double[] FvTarget;
                    string sourceFileT = modelPath + s + ".fv_B_" + B;
                    if (File.Exists(sourceFileT))
                    {
                        Console.WriteLine("Obtaining feature vector from file");
                        FvTarget = FileParser.parseFeatureVector(modelPath, s, B);
                    }
                    else
                    {
                        Console.WriteLine("Calculating feature vector");
                        FvTarget = extractFeatureVector(ref tarModel, ref directSorted);
                        FileWriter.WriteFeatureFile(FvTarget, modelPath, s);
                    }
                    Console.WriteLine("Comparing " + f + " AND " + s);
                    
                    distances[total] = GetDistance(ref FvCurrent, ref FvTarget);
                    total++;
                }


                FileWriter.WriteDistanceFile(distances, modelPath + f, files);
            }
            DateTime stopTimeTotal = DateTime.Now;
            TimeSpan totalDuration = stopTimeTotal - startTimeTotal;
            Console.WriteLine("Total write duration: ", totalDuration.TotalHours + " " + totalDuration.TotalMinutes + " " + totalDuration.TotalSeconds);

            Console.WriteLine("Done");

            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);
            Application.Run(new Form1());
        }

        //Get the distance between two vectors
        private static double GetDistance(ref double[] FvCurrent, ref double[] FvTarget)
        {

            //the paper suggests a L1 metric? page 137, then P = 1.
            double result = 0;
            for (int i = 0; i < FvCurrent.Length; i++)
            {
                result += System.Math.Abs(FvCurrent[i] - FvTarget[i]);
            }

            return result;
        }

        //creates a sorted array (in phi and theta) of direction vectors
        static unitVector[] getDirections()
		{
			unitVector[] direct = new unitVector[4 * B * B];
			int total = 0;

            double tempTheta = (System.Math.PI * 2) / (4 * B), tempTheta2 = System.Math.PI/(4*B);
            double tempPhi = System.Math.PI / B;

			for (int a = 0; a < B2; a++)
				for (int b = 0; b < B2; b++)
				{
					double theta, phi;
                    theta = a * tempTheta + tempTheta2;
					phi = b * tempPhi; //equals 2bPI/2B
					direct[total] = new unitVector(phi, theta, a, b,
													new Vec3D(System.Math.Cos(phi) * System.Math.Sin(theta),
																System.Math.Sin(phi) * System.Math.Sin(theta),
																System.Math.Cos(theta)));
					total++;
				}
			return direct;
		}

		//do all the important stuff :-)
		static double[] extractFeatureVector(ref Model3D model, ref unitVector[] ds)
		{
            //bruteforce
           // double[] extends2 = getExtendBruteForce(ref model, ref ds);
            
            //optimized
            double[] extends = getExtend(ref model, ref ds);

            //because of a error with the calculation of phi because of Atan2, the 0 values are recalculated
            for (int i = 0; i < extends.Length; i++)
            {
                if (extends[i] == 0.0)
                {
                    int b = i % B2;
                    int a = (i / B2);
                    extends[i] = getExtendOneVector(ref model, ref ds[B2 * b + a]);
                }
            }
            
            double[] feature = SFFT(extends, ds);
            return feature;
		}

        private static double getExtendOneVector(ref Model3D model, ref unitVector direction)
        {
            double tempExtend = 0.0;
            for (int i = 0; i < model.faces.Length; i++)
            {
                Face f = model.faces[i];
                double extend = RayTriangleIntersection2(ref direction, ref model, ref f);
                if (extend > tempExtend)
                    tempExtend = extend;
            }

            return tempExtend;
        }

        private static double[] switchOrder(double[] extends)
        {
            double[] ex = new double[extends.Length];

            for(int a=0; a < B2; a++)
                for (int b = 0; b < B2; b++)
                {
                    ex[a*B2+b] = extends[ b*B2+a ];
                }

            return ex;
        }

        static double[] getExtend(ref Model3D model, ref unitVector[] ds)
        {
            double[] extends = new double[4*B*B];

            //initialize all values to 0
            for(int q = 0; q < extends.Length; q++)
                extends[q] = 0.0;

            int counter = 0;

            for(int k=0; k < model.faces.Length; k++)
            {
                Face f = model.faces[k];
                double[] phiTheta = getMinMaxPhiTheta(ref f, ref model.vertices,ref model);
                double minPhi = phiTheta[0];
                double minTheta = phiTheta[1];
                double maxPhi = phiTheta[2];
                double maxTheta = phiTheta[3];

                double searchPhi;
                if (minPhi < 0)
                    searchPhi = minPhi + (2 * System.Math.PI);
                else
                    searchPhi = minPhi;

                int minPhiIndex = BinSearchPhi(ref ds, searchPhi, 0, ds.Length - 1, false);
                while (minPhiIndex > 0 && ds[minPhiIndex - 1].phi == ds[minPhiIndex].phi)
                    minPhiIndex--;

                if (maxPhi < 0)
                    searchPhi = maxPhi + 2 * System.Math.PI;
                else
                    searchPhi = maxPhi;

                int maxPhiIndex = BinSearchPhi(ref ds, maxPhi, 0, ds.Length - 1, true);
                while (maxPhiIndex < ds.Length-1 && ds[maxPhiIndex + 1].phi == ds[maxPhiIndex].phi)
                    maxPhiIndex++;

                //min en max van de theta zijn onveranderlijk over de subset in phi. Dit vanwege 
                //de sortering over eerst phi en vervolgens theta, waardoor de opeenvolging altijd hetzelfde is.
                int minThetaIndex = BinSearchTheta(ref ds, minTheta, minPhiIndex, minPhiIndex + B2, false);
                minThetaIndex -= minPhiIndex; // maakt de index relatief ten opzicht van de phiIndex

                int maxThetaIndex = BinSearchTheta(ref ds, maxTheta, minPhiIndex, minPhiIndex + B2, true);
                maxThetaIndex -= minPhiIndex;

                int thetaRange = maxThetaIndex - minThetaIndex;

                int range = System.Math.Abs(maxPhiIndex - minPhiIndex) / B2;

                for (int i = 0; i <= range; i++)
                {
                    for (int j = 0; j <= thetaRange; j++)
                    {
                        int index = (minPhiIndex + i * B2 + minThetaIndex + j) % ds.Length;
                        //Modulo vanwege dat als maxPhi groter is dan PI het eigenlijk negatief is, daarom moet er over het einde van de array worden doorgelopen.
                        unitVector temp = ds[index]; 
                        double extend = RayTriangleIntersection2(ref temp, ref model, ref f);
                        counter++;
                        if (extend == 0.0)
                        {
                           // System.Console.WriteLine("Did not hit. a: {0} b: {1}",temp.a, temp.b);
                            continue;
                        }
                        int temporary = (temp.a * 2 * B) + temp.b;
                        if(extends[temporary] < extend)
                            extends[temporary] = extend;

/*                        Console.WriteLine("Extends {0}: \textend: {1},\tphi: {2},\ttheta: {3},\t" +
                                          "minPhi: {4},\tmaxPhi: {5},\tminTheta: {6},\tmaxTheta: {7},\t" +
                                          "minPhiIndexValue: {8},\tmaxPhiIndexValue: {9},\tminThetaIndexValue: {10},\tmaxThetaIndexValue: {11}",
                                          temporary, extend, temp.phi, temp.theta, minPhi, maxPhi,
                                          minTheta, maxTheta, ds[minPhiIndex + minThetaIndex].phi, ds[minPhiIndex + (range - 1) * B2 + minThetaIndex].phi,
                                          ds[minPhiIndex + minThetaIndex].theta, ds[minPhiIndex + (range - 1) * B2 + minThetaIndex].theta);
  */                }
                }
            }

            System.Console.WriteLine(counter);
            
            return extends;
        }

        static double[] getExtendBruteForce(ref Model3D model, ref unitVector[] ds)
        {
            double[] extends = new double[4 * B * B];

            //initialize all values to 0
            for (int q = 0; q < extends.Length; q++)
                extends[q] = 0.0;

            int counter = 0;

            for (int k = 0; k < model.faces.Length; k++)
            {
                Face f = model.faces[k];

                for (int i = 0; i < ds.Length; i++)
                {
                    unitVector temp = ds[i];
                    double extend = RayTriangleIntersection2(ref temp, ref model, ref f);
                    counter++;
                    if (double.IsNaN(extend))
                    {
                       // System.Console.WriteLine("Did not hit. a: {0} b: {1}", temp.a, temp.b);
                        continue;
                    }
                    int temporary = (temp.a * 2 * B) + temp.b;
                    if (extends[temporary] < extend)
                        extends[temporary] = extend;
    /*                Console.WriteLine("Extends {0}:\textend: {1},\tphi: {2},\ttheta: {3},\t" +
                                          "minPhi: {4},\tmaxPhi: {5},\tminTheta: {6},\tmaxTheta: {7},\t" +
                                          "minPhiIndexValue: {8},\tmaxPhiIndexValue: {9},\tminThetaIndexValue: {10},\tmaxThetaIndexValue: {11}",
                                          temporary, extend, temp.phi, temp.theta, 0, System.Math.PI*2,
                                          0, System.Math.PI, 0, System.Math.PI * 2,
                                          0, System.Math.PI);
      */        }
            }

            return extends;
        }

        static int BinSearchTheta(ref unitVector[] ds, double value, int low, int high, bool max)
        {
            if (high == low)
                return low;

            if (high < low)
                if (max)
                    return low;
                else
                    return high;

            int mid = low + (int)((high - low) / 2);
            if (ds[mid].theta > value)
                return BinSearchTheta(ref ds, value, low, mid - 1, max);
            else if (ds[mid].theta < value)
                return BinSearchTheta(ref ds, value, mid + 1, high, max);
            else
                return mid; // found
        }

        static int BinSearchPhi(ref unitVector[] ds, double value, int low, int high, bool max)
        {
            if (high == low)
                return low;
            int mid = low + (int)((high - low) / 2);
            if (ds[mid].phi > value)
               return BinSearchPhi(ref ds, value, low, mid - 1, max);
            else if (ds[mid].phi < value)
               return BinSearchPhi(ref ds, value, mid + 1, high, max);
            else
               return mid; // found

        }

        /**
		 * Methode die het punt bepaald waarop een ray, gerepresenteerd door a en b in een array, in een face f. 
		 * De vertices van de Face f zitten in de model waar f in zit.
		 */
        static double RayTriangleIntersection2(ref unitVector direction, ref Model3D model, ref Face f)
        {
            double r = 0.0;
            double d = 0.0;
            Vec3D[] v = model.vertices;
            Vec3D u = direction.vector;

            Vec3D A = v[f.vertices[0]];
            Vec3D B = v[f.vertices[1]];
            Vec3D C = v[f.vertices[2]];
            Vec3D CA = A.subtract(C);
            Vec3D CB = B.subtract(C);

            Vec3D p = u.cross(CB);
            double det = CA.dot(p);
            if (det > double.Epsilon)
            {
                double AlphaS = C.times(-1.0).dot(p);
                if (AlphaS >= 0 && AlphaS <= det)
                {
                    Vec3D q = CA.cross(C);
                    double BetaS = u.dot(q);
                    if (BetaS >= 0 && AlphaS + BetaS <= det)
                        d = (CB.dot(q)) / det;
                }
            }
            else if(det < 0-Double.Epsilon)
            {
                double AlphaS = C.times(-1).dot(p);
                if (AlphaS <= 0 && AlphaS >= det)
                {
                    Vec3D q = CA.cross(C);
                    double BetaS = u.dot(q);
                    if (BetaS <= 0 && AlphaS + BetaS >= det)
                        d = (CB.dot(q)) / det;
                }
            }

            if (d > r)
                r = d;
            //the ray will not intersect with the triangle
            return r;

        }

		/**
		 * Methode die het punt bepaald waarop een ray, gerepresenteerd door a en b in een array, in een face f. 
		 * De vertices van de Face f zitten in de model waar f in zit.
		 */
		static double RayTriangleIntersection(ref unitVector direction, ref Model3D model, ref Face f)
		{
			double r = 0.0;
			Vec3D[] v = model.vertices;

			Vec3D h = f.normal;// v[f.vertices[1]].subtract(v[f.vertices[0]]).cross(v[f.vertices[2]].subtract(v[f.vertices[0]]));
			if (h.length() > double.Epsilon)
			{
				Vec3D normal = f.normal.normalize();
				//normalize the normal, then dot it with the direction vector
				double NU = System.Math.Abs(normal.dot(direction.vector));
				if (NU > double.Epsilon)
				{
					double d = (normal.dot(v[f.vertices[0]])) / NU;
					if (d > r)
					{
						Vec3D P = direction.vector.times(d);
						//is P on the triangle?
						if ((v[f.vertices[0]].subtract(P).cross(v[f.vertices[1]].subtract(P))).dot(f.normal) >= 0 &&
							(v[f.vertices[1]].subtract(P).cross(v[f.vertices[2]].subtract(P))).dot(f.normal) >= 0 &&
							(v[f.vertices[2]].subtract(P).cross(v[f.vertices[0]].subtract(P))).dot(f.normal) >= 0)
							return d;
					}

				}
			}

			//the ray will not intersect with the triangle
			return double.NaN;

		}

        static double zRound(double z)
        {
            while (z < 0.0) z += 2 * System.Math.PI;
            while (z > 2 * System.Math.PI) z -= 2 * System.Math.PI;

            return z;
        }

		static double[] getMinMaxPhiTheta(ref Face f, ref Vec3D[] v, ref Model3D model)
		{
            /*	int minPhi      = phiTheta[0];
				int minTheta    = phiTheta[1];
				int maxPhi      = phiTheta[2];
				int maxTheta    = phiTheta[3]; */
            double[] minmax = new double[4];
            minmax[0] = double.MaxValue;
            minmax[1] = double.MaxValue;
            minmax[2] = double.MinValue;
            minmax[3] = double.MinValue;

            double  minX = double.MaxValue, maxX = double.MinValue, minY = double.MaxValue,
                    maxY = double.MinValue, minZ = double.MaxValue, maxZ = double.MinValue, meanX = 0, meanY = 0, meanZ = 0;

            for(int k=0; k < f.numberOfVertices; k++)
            {
                Vec3D vec = v[f.vertices[k]];
                updateMinMaxMeanXYZ(ref minX,ref maxX,ref meanX,vec.x,ref minY,ref maxY,ref meanY,vec.y,ref minZ,ref maxZ,ref meanZ,vec.z,k);
                double phi, theta;
                phi = System.Math.Atan2(vec.y, vec.x);
                theta = System.Math.Acos((vec.z) / (vec.length()));

                if (minmax[0] > phi)
                    minmax[0] = phi;
                if (minmax[1] > theta)
                    minmax[1] = theta;
                if (minmax[2] < phi)
                    minmax[2] = phi;
                if (minmax[3] < theta)
                    minmax[3] = theta;
            }

            /*
            //Because of the possibility of an polygon lying on all equal or close theta, the correct theta may not be found with only the points.
            //Trying to solve this here.
            Vec3D dotPoint = new Vec3D(double.NaN,double.NaN,double.NaN);
            double dotDist = double.MaxValue;

            for (int k = 0; k < f2.numberOfVertices; k++)
            {
                Vec3D a, b, edge;
                a = v[f.vertices[k]].copy().normalize();
                b = v[f.vertices[(k + 1) % f.numberOfVertices]].copy().normalize();
                edge = a.subtract(b).normalize();

                double dotO = b.dot(edge);

                if (dotO < 0.0)
                {
                    Vec3D tempPoint = b.addition(edge.times(dotO));
                    double tempDist = tempPoint.squaredLength();
                    if (dotDist > tempDist)
                    {
                        dotDist = tempDist;
                        dotPoint = tempPoint;
                    }
                }
            }

            minmax[3] = System.Math.Acos((dotPoint.z) / (dotPoint.length()));
            */

            //Check om te bepalen of de polygon in het xy vlak (z is up),
            //zo ja, check of erboven of eronder is en pas theta aan om volledig door te lopen.
            if (minX < 0 && maxX > 0 && minY < 0 && maxY > 0)
                if (meanZ > 0)
                    minmax[1] = 0;
                else
                    minmax[3] = System.Math.PI;

            return minmax;
		}



        private static void updateMinMaxMeanXYZ(ref double minX, ref double maxX, ref double meanX, double x,
                                                ref double minY, ref double maxY, ref double meanY, double y,
                                                ref double minZ, ref double maxZ, ref double meanZ, double z, int k)
        {
            if (x < minX)
                minX = x;
            if (x > maxX)
                maxX = x;
            meanX = (meanX*k+x)/(k+1); 

            if (y < minY)
                minY = y;
            if (y > maxY)
                maxY = y;
            meanY = (meanY * k + y) / (k + 1); 

            if (z < minZ)
                minZ = z;
            if (z > maxZ)
                maxZ = z;
            meanZ = (meanZ * k + z) / (k + 1); 
        }

        static public double[] SFFT(double[] data, unitVector[] ds)
        {
            //B^2 complex coefficients Flm in the range 0 (<= |m| <= l <= B)
            //extract the descriptor for k! the paper says 29 is good, or use the rotation invariant feature vector

            double[] rcoeffs = new double[B * B];   //(double*)malloc(sizeof(double) * (bw * bw));
            double[] icoeffs = new double[B * B];   //(double*)malloc(sizeof(double) * (bw * bw));
            double[] rdata = data;                  //(double*)malloc(sizeof(double) * (size * size));
            double[] idata = new double[4 * B * B];     //(double*)malloc(sizeof(double) * (size * size));
            double[] rdata2 = new double[4 * B * B];                //(double*)malloc(sizeof(double) * (size * size));
            double[] idata2 = new double[4 * B * B];     //(double*)malloc(sizeof(double) * (size * size));
            //idata = 0.0
            for (int i = 0; i < idata.Length; i++)
            {
                idata[i] = 0.0;
            }

            Console.WriteLine("Starting SFFT");
            DateTime startTime1 = DateTime.Now;

            try
            {
                calcCoeff(B, B, rdata, idata, rcoeffs, icoeffs);

                //if you wanna compute the real model again
                inverseTransform(B, B, rcoeffs, icoeffs, rdata2, idata2);
            }
            catch (DllNotFoundException e)
            {
                Console.WriteLine(e.StackTrace + "\n");
                Console.WriteLine("DLL not found.");
                Application.Exit();
            }
            //if you wanna compute the real model again
            
            Vec3D[] points = new Vec3D[4*B*B];
            for (int i = 0; i < rdata2.Length; i++)
            {
                int b = i % B2;
                int a = (i / B2);
                points[i] = ds[B2*b+a].vector.times(rdata2[i]);
            }
            Random test = new Random();
            FileWriter.WritePointCloudObj(points, "pointcloud" + test.Next(100));
            

            DateTime stopTime1 = DateTime.Now;
            TimeSpan duration1 = stopTime1 - startTime1;
            Console.WriteLine("SFFT duration: {0}", duration1.TotalSeconds);
            Console.WriteLine("Done SFFT");

            double[] fv = new double[B];
            /*
            int k = 29;
            double[] fv2 = new double[k * (k + 1) / 2 + 1];
            int total = 0;
            for (int l = 0; l < k; l++)
            {
                for (int m = 0; m < l; m++)
                {
                    int index = getCoeffIndex(m, l);
                    double mag = System.Math.Sqrt(rcoeffs[index] * rcoeffs[index] + icoeffs[index] * icoeffs[index]);
                    if(m == 0)
                    {
                        fv2[total] = mag/2;
                        total++;
                    }
                    else
                    {
                        fv2[total] = mag;
                        total++;
                    }
                   
                }
            }

            return fv2; */

            for (int l = 0; l < B; l++)
            {
                double result = 0;
                for (int m = 0 - l; m <= l; m++)
                {
                    int index = getCoeffIndex(m, l);
                    result += rcoeffs[index] * rcoeffs[index] + icoeffs[index] * icoeffs[index];
                }
                fv[l] = System.Math.Sqrt(result);
            }

            return fv;
        }

        //get the right index from the coeff arrays, leeched from the C function
        private static int getCoeffIndex(int m, int l)
        {
            int bigL = B - 1;

            if (m >= 0)
                return (m * (bigL + 1) - ((m * (m - 1)) / 2) + (l - m));
            else
                return (((bigL * (bigL + 3)) / 2) + 1 +
                    ((bigL + m) * (bigL + m + 1) / 2) + (l - System.Math.Abs(m)));
        }

	}
}