import math
import scipy.optimize as opt

def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa

    Return: העומס P בניוטון (float)
    """
    
    # מומנט האינרציה (I = A * r^2)
    I = A * (r ** 2)
    
    # חישוב עומס אוילר המשמש כחסם עליון לחיפוש שלנו
    euler_limit = (math.pi ** 2 * E * I) / (L ** 2)
    
    # חישוב קבועים מראש כדי למנוע חישובים מיותרים בתוך הלולאה
    eccentricity_ratio = (e * c) / (r ** 2)
    
    # הפונקציה שאותה אנחנו רוצים לאפס (הפרש המאמצים)
    def equation_to_solve(current_p):
        if current_p <= 0:
            return -sigma_allow
            
        # הארגומנט של פונקציית הסקאנט
        angle_rad = (L / (2 * r)) * math.sqrt(current_p / (E * A))
        secant_val = 1.0 / math.cos(angle_rad)
        
        # חישוב המאמץ הנוכחי
        current_stress = (current_p / A) * (1.0 + eccentricity_ratio * secant_val)
        
        # מר
