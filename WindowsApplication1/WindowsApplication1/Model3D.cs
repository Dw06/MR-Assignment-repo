/*
 * Multimedia Retrieval Practical Assignment: Ray-Based Approach with Spherical Harmonics
 * By Milan de Graaf, 3117308 and David Weterings, 3117480
 * 
 **/
namespace WindowsApplication1
{
    //Class that represents a 3D-model
    class Model3D
    {
        public Vec3D[] vertices;
        public Face[] faces;
        public Vec3D COM;

        public int numberOfVertices;
        public int numberOfFaces;

        public string source;
        public double totalArea;

        public Model3D(Vec3D[] v, Face[] f, string s)
        {
            vertices = v;
            faces = f;
            numberOfVertices = v.Length;
            numberOfFaces = f.Length;
            source = s;
            totalArea = 0.0;
            foreach (Face fac in f)
                totalArea += fac.area;
            this.TransToCOM();
            //this.normalizeScale();
        }

        //get CoM, if it doesn't exist: set it!
        public Vec3D getCOM()
        {
            if (COM != null)
                return COM;

            double x = 0, y = 0, z = 0;
            foreach (Vec3D v in vertices)
            {
                x += v.x;
                y += v.y;
                z += v.z;
            }
            COM = new Vec3D(x / vertices.Length, y / vertices.Length, z / vertices.Length);
            return COM;
        }

        //translate every vertex relative to CoM
        public void TransToCOM()
        {
            Vec3D mean = getCOM();

            if (COM.Equals(new Vec3D(0, 0, 0)))
            {
                return;
            }

            for (int i = 0; i < vertices.Length; i++)
            {
                Vec3D v = vertices[i];
                v.sub(mean);
            }
            COM = new Vec3D(0, 0, 0);
        }

        public void normalizeScale2()
        {
            Vec3D CoM = getCOM();

            double distance = 0;

            foreach (Vec3D v in vertices)
            {
                double x = v.x - COM.x;
                double y = v.y - COM.y;
                double z = v.z - COM.z;

                distance += getDistanceOfPoint(x, y, z);
            }

            distance /= vertices.Length;
        }

        private double getDistanceBetweenPoints(Vec3D v1, Vec3D v2)
        {
            return System.Math.Sqrt(System.Math.Pow(v1.x - v2.x, 2) + System.Math.Pow(v1.y - v2.y, 2) + System.Math.Pow(v1.z - v2.z, 2));
        }

        private double getDistanceBetweenPoints(double x1, double y1, double z1, double x2, double y2, double z2)
        {
            return System.Math.Sqrt(System.Math.Pow(x1 - x2, 2) + System.Math.Pow(y1 - y2, 2) + System.Math.Pow(z1 - z2, 2));
        }

        private double getDistanceOfPoint(double x, double y, double z)
        {
            return System.Math.Sqrt(x * x + y * y + z * z);
        }

        private double getDistanceOfPoint(Vec3D v)
        {
            return System.Math.Sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        }

        private void normalizeScale()
        {
            double Sx = 0.0, Sy = 0.0, Sz = 0.0;
            for (int i = 0; i < numberOfFaces; i++)
            {
                Face f = faces[i];
                double faceArea = f.area;
                Sx += calcM(ref f.vertices, 'x') + faceArea;
                Sy += calcM(ref f.vertices, 'y') + faceArea;
                Sz += calcM(ref f.vertices, 'z') + faceArea;
            }

            Sx *= (1 / (3 * totalArea));
            Sy *= (1 / (3 * totalArea));
            Sz *= (1 / (3 * totalArea));
            double invScaleFactor = System.Math.Pow(System.Math.Sqrt((Sx*Sx+Sy*Sy + Sz*Sz)/3),-1.0);
            foreach (Vec3D v in vertices)
                v.tim(invScaleFactor);

        }

        private double calcM(ref int[] f, char test)
        {
            double Xa, Xb, Xc;
            if(test == 'x')
            {
                Xa = vertices[f[0]].x; 
                Xb = vertices[f[1]].x;
                Xc = vertices[f[2]].x;
            }
            else if(test == 'y')
            {
                Xa = vertices[f[0]].y; 
                Xb = vertices[f[1]].y;
                Xc = vertices[f[2]].y;
            }
            else
            {
                Xa = vertices[f[0]].z; 
                Xb = vertices[f[1]].z;
                Xc = vertices[f[2]].z;
            }

            if (Xa * Xb >= 0 &&
                Xb * Xc >= 0 &&
                Xa * Xc >= 0)
                return System.Math.Abs(Xa + Xb + Xc);
            else if ((Xa * Xb <= 0 && Xa * Xc <= 0) || (Xb * Xc >= 0))
                return System.Math.Abs(Xa + Xb + Xc) - 2 * (Xa * Xa * Xa) / ((Xb - Xa) * (Xc - Xa));
            else
                return 0.0;
        }

        /*
        private void CPCA()
        {
            //invariant(pointset) = s^-1 * F * A * (v - Mi)
            //Si = area van face i, S = total area
            //center
            Vec3D Mi = new Vec3D(0, 0, 0);
            for (int m = 0; m < numberOfFaces; m++)
            {
                Mi.add(faces[m].CoG.times(faces[m].area));
            }
            Mi.divide(this.totalArea);
            double[,] Ci = new double[numberOfVertices,numberOfVertices];




        } */
    }
}
