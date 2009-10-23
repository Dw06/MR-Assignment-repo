/*
 * Multimedia Retrieval Practical Assignment: Ray-Based Approach with Spherical Harmonics
 * By Milan de Graaf, 3117308 and David Weterings, 3117480
 * 
 **/
using System;
using System.Collections;

namespace WindowsApplication1
{
    //Contains a few generic classes and small datastructures needed for this project.


    public class Vec3D
    {
        public double x, y, z;

        public Vec3D()
        {
            x = y = z = 0.0;
        }

        public Vec3D(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        public Vec3D copy()
        {
            return new Vec3D(x, y, z);
        }

        public void add(Vec3D v)
        {
            this.x += v.x;
            this.y += v.y;
            this.z += v.z;
        }

        public void tim(double d)
        {
            this.x *= d;
            this.y *= d;
            this.z *= d;
        }

        public Vec3D addition(Vec3D v)
        {
            return new Vec3D(this.x + v.x, this.y + v.y, this.z + v.z);
        }

        public void div(double d)
        {
            this.x /= d;
            this.y /= d;
            this.z /= d;
        }

        public Vec3D divide(double d)
        {
            return new Vec3D(this.x / d, this.y / d, this.z / d);
        }

        public void sub(Vec3D v)
        {
            this.x -= v.x;
            this.y -= v.y;
            this.z -= v.z;
        }

        public Vec3D times(double f)
        {
            return new Vec3D(f * x, f * y, f * z);
        }

        public Vec3D subtract(Vec3D v)
        {
            return new Vec3D(x - v.x, y - v.y, z - v.z);
        }

        public Vec3D subtract(double v)
        {
            return new Vec3D(x - v, y - v, z - v);
        }

        public double squaredLength()
        {
            return x * x + y * y + z * z;
        }

        public double length()
        {
            return System.Math.Sqrt(squaredLength());
        }

        public Vec3D normalize()
        {
            double l = this.length();
            if (l != 0.0)
            {
                return new Vec3D(x /= l, y /= l, z /= l);
            }
            return new Vec3D();
        }

        public double dot(Vec3D v)
        {
            return x * v.x + y * v.y + z * v.z;
        }

        public Vec3D cross(Vec3D v)
        {
            return new Vec3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
        }

        public bool Equals(Vec3D v)
        {
            // If parameter is null return false:
            if ((object)v == null)
            {
                return false;
            }

            // Return true if the fields match:
            return (x == v.x) && (y == v.y) && (z == v.z);
        }
    }

    public class Face
    {
        public int numberOfVertices;
        public int[] vertices;
        public Vec3D normal;
        public double area;

        public Face(int[] v, Vec3D normal, double area)
        {
            vertices = v;
            numberOfVertices = v.Length;
            this.normal = normal;
            this.area = area;
        }
    }

    //datastructure used for the directions on the unit sphere
    public class unitVector : IComparable
    {
        public double phi, theta;
        public Vec3D vector;
        public int a, b;

        public unitVector(double phi, double theta, int a, int b, Vec3D v)
        {
            this.phi = phi;
            this.theta = theta;
            this.a = a;
            this.b = b;
            vector = v;
        }

        public int CompareTo(object obj)
        {
            unitVector u = (unitVector)obj;

            if (this.phi < u.phi || (this.phi == u.phi && this.theta < u.theta))
            {
                //this < u
                return -1;
            }
            else if (this.phi == u.phi && this.theta == u.theta)
            {
                //equal
                return 0;
            }
            else
                return 1;
        }

    }

    //test class
    public class TestTuple : IComparable
    {
        int number;
        double extends, phi, theta;

        public TestTuple(int number, double extends, double phi, double theta)
        {
            this.number = number;
            this.extends = extends;
            this.phi = phi;
            this.theta = theta;
        }

        public int CompareTo(object obj)
        {
            unitVector u = (unitVector)obj;

            if (this.phi < u.phi || (this.phi == u.phi && this.theta < u.theta))
            {
                //this < u
                return -1;
            }
            else if (this.phi == u.phi && this.theta == u.theta)
            {
                //equal
                return 0;
            }
            else
                return 1;
        }
    }
}
