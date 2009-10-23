/*
 * Multimedia Retrieval Practical Assignment: Ray-Based Approach with Spherical Harmonics
 * By Milan de Graaf, 3117308 and David Weterings, 3117480
 * 
 **/
using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Text;

namespace WindowsApplication1
{
    //FileParser is responsible for opening/parsing all files for this project
    static class FileParser
    {

        public static Dictionary<string, ArrayList> parseBaseFile(string modelPath, string sFilename)
        {
            StreamReader sr = new StreamReader(modelPath + sFilename);
            string line;
            Dictionary<string, ArrayList> categories = new Dictionary<string, ArrayList>();

            while ((line = sr.ReadLine()) != null)
            {
                string[] items = line.Split('/');
                string[] extension = items[items.Length - 1].Split('.');
                Model3D m3d;
                if (extension.Length > 1)
                {
                    if (extension[extension.Length - 1] == "obj")
                    {
                        m3d = parseObjModel(modelPath + line, line);
                    }
                    else //if extension unknown
                        throw new Exception("This is a unknown data format.");
                }
                else
                    m3d = parseSofModel(modelPath + line, line);

                if (categories.ContainsKey(items[0]))
                {
                    categories[items[0]].Add(m3d);
                }
                else
                {
                    ArrayList m3dArr = new ArrayList();
                    m3dArr.Add(m3d);
                    categories.Add(items[0], m3dArr);
                }
            }
            sr.Close();
            return categories;
        }

        //Parse the basenames file, but in a simple way, line for line
        public static ArrayList parseBaseFileSimple(string modelPath, string sFilename)
        {
            StreamReader sr = new StreamReader(modelPath + sFilename);
            string line;
            ArrayList files = new ArrayList();

            while ((line = sr.ReadLine()) != null)
            {
                files.Add(line);
            }
            sr.Close();
            return files;
        }

        //Parse a .obj file into a Model3D
        private static Model3D parseObjModel(string path, string sourcepath)
        {
            StreamReader sr = new StreamReader(path);
            string line = sr.ReadLine();
            ArrayList vertices = new ArrayList();
            ArrayList faces = new ArrayList();
            
            while (line != null)
            {
                if (line.Length != 0)
                {
                    //de obj opdracht 'vn' en 'vt' kijken we niet naar, alleen 'v' is van belang.
                    if (line[0] == 'v')
                    {
                        if (line[1] == ' ')
                        {
                            vertices.Add(parseObjVertex(line.Substring(3)));
                        }
                    }

                    //allen de opdracht 'f' is van belang, weet niet of er nog andere opdrachten met f beginnen, maar voor robuustheid op deze manier afgehandelt.
                    //f komt niet voor voordat alle vertices zijn afgehandelt.
                    else if (line[0] == 'f')
                    {
                        if (line[1] == ' ')
                        {
                            faces.AddRange(parseObjFace(line.Substring(2), ref vertices));
                        }
                    }
                }


                line = sr.ReadLine();
            }

            sr.Close();

            Vec3D[] v = (Vec3D[])vertices.ToArray(typeof(Vec3D));
            Face[] f = (Face[])faces.ToArray(typeof(Face));

            return new Model3D(v, f, sourcepath);
        }

        //parse the face for an .obj file
        private static ArrayList parseObjFace(string line, ref ArrayList vertices) //Arraylist van Face is de return
        {
            ArrayList faces = new ArrayList();
            ArrayList verts = new ArrayList();
            line = line.Trim();
            string[] vertsInts = line.Split(' ');
            //NOTE THE -1 FOR CORRECT INDICES
            foreach (string s in vertsInts)
                verts.Add(System.Convert.ToInt32(s)-1);

            int[] v = (int[])verts.ToArray(typeof(int));

            Vec3D[] ver = (Vec3D[])vertices.ToArray(typeof(Vec3D));

            if (v.Length > 3)
                faces = Triangulate(ref ver, new Face(v, new Vec3D(0, 0, 0), 0));
            else
            {
                Vec3D a, b, c;
                a = (Vec3D)vertices[v[0]];
                b = (Vec3D)vertices[v[1]];
                c = (Vec3D)vertices[v[2]];

                Vec3D normal = (b.subtract(a)).cross(c.subtract(a));
                double area = 0.5 * System.Math.Abs((c.subtract(a).cross(b.subtract(a)).length()));
                faces.Add(new Face(v, normal, area));
            }

            return faces;
        }

        //parse the vertex for an .obj file
        private static Vec3D parseObjVertex(string line)
        {
            line = line.Trim();
            //assumption that there is only 1 whitespace inbetween the numbers.
            string[] positions = line.Split(' ');
            double x = System.Convert.ToDouble((positions[0]));
            double y = System.Convert.ToDouble((positions[1]));
            double z = System.Convert.ToDouble((positions[2]));
            return new Vec3D(x, y, z);
        }

        //parse een SOF file
        public static Model3D parseSofModel(string path, string sourcepath)
        {
            StreamReader sr = new StreamReader(path);

            string line;
            int vertices, faces;
            vertices = System.Convert.ToInt32(sr.ReadLine());
            faces = System.Convert.ToInt32(sr.ReadLine());
            Vec3D[] v = new Vec3D[vertices];
            ArrayList f = new ArrayList();

            for (int i = 0; i < vertices; i++) //possibly more than one whitespace inbetween, use Scanner?
            {
                line = sr.ReadLine();
                string[] items = line.Split(' ');
                try
                {
                    v[i] = new Vec3D(System.Convert.ToDouble(items[0]),
                                    System.Convert.ToDouble(items[1]),
                                    System.Convert.ToDouble(items[2]));
                }
                catch(FormatException)
                {
                    return null;
                }
            }

            for (int i = 0; i < faces; i++) //possibly more than one whitespace inbetween, use Scanner?
            {
                line = sr.ReadLine();
                string[] items = line.Split(' ');
                int aantal = System.Convert.ToInt32(items[0]);
                int[] vIndex = new int[aantal];
                for (int z = 1; z <= aantal; z++)
                {
                    vIndex[z - 1] = System.Convert.ToInt32(items[z]);
                }

                //Assumption of a polygon lying on a plane, would be hard otherwise, an average over the normals on the surface spline.
                Vec3D a, b, c;
                a = v[vIndex[0]];
                b = v[vIndex[1]];
                c = v[vIndex[2]];

                Vec3D normal = (b.subtract(a)).cross(c.subtract(a));

                ArrayList temp = Triangulate(ref v, new Face(vIndex, normal,0));
                foreach (Face face in temp)
                    f.Add(face);
            }
            sr.Close();

            Face[] f2 = (Face[])f.ToArray(typeof(Face));

            return new Model3D(v, f2, sourcepath);
        }

        static ArrayList Triangulate(ref Vec3D[] v, Face f)
        {
            //verdeel Face f in meerdere kleine faces van 3 vertices
            //Er wordt vanuit gegaan dat een face meer dan drie vertices heeft als deze in deze methode komt.
            //en geef deze in een arraylist terug, MOGELIJK in een array.

            //BEGIN Naive/Simple implementation, ONLY works on convex polygons.
            ArrayList q = new ArrayList();


            for (int i = 1; i < f.numberOfVertices - 1; i++)
            {
                int[] verts = new int[3];
                verts[0] = f.vertices[0];
                verts[1] = f.vertices[i];
                verts[2] = f.vertices[i + 1];

                Vec3D a, b, c;
                a = v[verts[0]];
                b = v[verts[1]];
                c = v[verts[2]];

                Vec3D norm = (b.subtract(a)).cross(c.subtract(a));
                double area = 0.5 * System.Math.Abs((c.subtract(a).cross(b.subtract(a)).length()));
                q.Add(new Face(verts, norm,area));
            }
            return q;

            //END Naive/Simple implementation
        }

        //parse a Model given a file, checks the extension
        internal static Model3D parseModel(string modelPath, string fileinfo)
        {
            string[] items = fileinfo.Split('/');
            string[] extension = items[items.Length - 1].Split('.');
            Model3D m3d;
            if (extension.Length > 1)
            {
                if (extension[extension.Length - 1] == "obj")
                {
                    m3d = parseObjModel(modelPath, fileinfo);
                }
                else //if extension unknown
                    throw new Exception("This is a unknown data format.");
            }
            else
                m3d = parseSofModel(modelPath, fileinfo);

            return m3d;
        }
        
        // Parse a feature vector file for a certain bandwidth
        internal static double[] parseFeatureVector(string modelPath, string f, int B)
        {
            StreamReader sr = new StreamReader(modelPath + f + ".fv_B_" + B);
            double[] fv = new double[B];

            for (int i = 0; i < B; i++)
                fv[i] = System.Convert.ToDouble(sr.ReadLine());

            return fv;
        }
    }
}