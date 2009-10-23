import jigl.image.GrayImage;

public class houghLine {
	double theta,r;
	int pixelvalue;
	
	houghLine(int t,int r, int p) {
		theta = t;
		this.r = r;
		pixelvalue = p;
	}

	
	public String toString() {
		return "HOUGHLINE: theta: " + theta + " r: " + r + " pxl: " + pixelvalue + "\n";
	}


	public void draw(GrayImage temp,int radius) {
		int height = temp.Y(); 
		int width = temp.X(); 
		double realTheta = Math.toRadians(theta-90);

		double sinT = Math.sin(realTheta); 
        double cosT = Math.cos(realTheta); 
        if (realTheta < Math.PI * 0.25 || realTheta > Math.PI * 0.75) { 
            //Verticale lijnen 
            for (int y = 0; y < height; y++) { 
                int x = (int) ((r + (y * sinT)) / cosT); 
                if (x < width && x >= 0) { 
                    temp.set(x, y, 0); 
                } 
            } 
        } else { 
            //Horizontaal
            for (int x = 0; x < width; x++) { 
                int y = (int) (r - (x * cosT) / sinT); 
                if (y < height && y >= 0) { 
                	temp.set(x, y, 0); 
                } 
            } 
        }
	}
	
	public void drawMax(GrayImage img) {
		for(int dt = (int) (theta-4); dt < theta+5; dt++) {
			for(int dr = (int) (r-5); dr < r+5; dr++) {
				if(dr>0 && dt>0) {
					img.set(dt, dr, 255);
				}
			}
		}
	}
}
