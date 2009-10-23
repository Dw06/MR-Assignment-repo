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
    //FileWriter is responsible for writing result data for this project
    static class FileWriter
    {
        /// <summary>
        /// Writes a distance file for a certain model
        /// </summary>
        /// <param name="distances"></param>
        /// <param name="targetDirectory"></param>
        /// <param name="files"></param>
        internal static void WriteDistanceFile(double[] distances, string targetDirectory, ArrayList files)
        {
            StreamWriter sw = new StreamWriter(targetDirectory + "._dist");
            int total = 0;
            foreach (string fileName in files)
            {
                sw.WriteLine(distances[total] + " " + fileName);
                total++;
            }
            sw.Close();
        }

        /// <summary>
        /// Writes a file for a model containing the feature vector
        /// </summary>
        /// <param name="fv"></param>
        /// <param name="modelPath"></param>
        /// <param name="f"></param>
        internal static void WriteFeatureFile(double[] fv, string modelPath, string f)
        {

            StreamWriter sw = new StreamWriter(modelPath + f + ".fv_B_" + fv.Length);
            foreach (double d in fv)
                sw.WriteLine(d);

            sw.Close();
        }
        /// <summary>
        /// Write a pointcloud for a vector array in txt format
        /// </summary>
        /// <param name="points"></param>
        /// <param name="filename"></param>
        internal static void WritePointCloudTxt(Vec3D[] points, string filename)
        {
            StreamWriter sw = new StreamWriter(filename + ".txt");
            foreach (Vec3D v in points)
            {
                sw.WriteLine(v.x + " " + v.y + " " + v.z);
            }
            sw.Close();
        }

        /// <summary>
        /// Write a pointcloud for a vector array in .obj format
        /// </summary>
        /// <param name="points"></param>
        /// <param name="filename"></param>
        internal static void WritePointCloudObj(Vec3D[] points, string filename)
        {
            StreamWriter sw = new StreamWriter(filename + ".obj");
            foreach (Vec3D v in points)
            {
                sw.WriteLine("v  " + v.x + " " + v.y + " " + v.z);
            }
            sw.Close();
        }
    }
}
