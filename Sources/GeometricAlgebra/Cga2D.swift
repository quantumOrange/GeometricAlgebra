import Foundation
import simd
// swift v0.002 Template for the preamble
public struct CGA2D  {
    var values:[Double]
    
    public static var zero:Self {
        CGA2D(values: Array(repeating: 0, count: 16))
    }
}

extension CGA2D {
    struct Vector {
        let  e1:Double
        let  e2:Double
        let  e:Double
        let  ebar:Double
        
        
        static func dot(_ a:Vector, _ b:Vector) -> Double {
            a.e1 * b.e1 +   a.e2 * b.e2 +   a.e * b.e -  a.ebar * b.ebar
        }
        
        static func * (_ a:Double, _ b:Vector) -> Vector {
            Vector(e1: a * b.e1, e2: a * b.e2, e: a * b.e, ebar: a * b.ebar)
        }
        
        var multivector:CGA2D {
            CGA2D(e1: e1, e2: e2, e: e, ebar: ebar)
        }
    }
    
    init( value:Double, index:Int) {
        var arr:[Double] = Array(repeating: 0.0, count: 16)
        arr[index] = value
        values = arr
    }
    
    init( e1 x:Double, e2 y:Double,e z:Double,ebar w:Double) {
        var arr:[Double] = Array(repeating: 0.0, count: 16)
        arr[1] = x
        arr[2] = y
        arr[3] = z
        arr[4] = w
        values = arr
    }
    
    // basis vectors are available
    public static let  e1 = CGA2D(value:1.0,index: 1)
    public static let  e2 = CGA2D(value:1.0,index: 2)
    public static let  e3 = CGA2D(value:1.0,index: 3)
    public static let  e4 = CGA2D(value:1.0,index: 4)
    public static let  e12 = CGA2D(value:1.0,index: 5)
    public static let  e13 = CGA2D(value:1.0,index: 6)
    public static let  e14 = CGA2D(value:1.0,index: 7)
    public static let  e23 = CGA2D(value:1.0,index: 8)
    public static let  e24 = CGA2D(value:1.0,index: 9)
    public static let  e34 = CGA2D(value:1.0,index: 10)
    public static let  e123 = CGA2D(value:1.0,index: 11)
    public static let  e124 = CGA2D(value:1.0,index: 12)
    public static let  e134 = CGA2D(value:1.0,index: 13)
    public static let  e234 = CGA2D(value:1.0,index: 14)
    public static let  e1234 = CGA2D(value:1.0,index: 15)
    
    public static let  e = CGA2D(value:1.0,index: 3)
    public static let  ebar = CGA2D(value:1.0,index: 4)
    public static let  n = CGA2D(value:1.0,index: 3) + CGA2D(value:1.0,index: 4) // e + ebar
    public static let  ñ = CGA2D(value:1.0,index: 3) - CGA2D(value:1.0,index: 4) // e - ebar
    
    public static func lift(point q:SIMD2<Double>) -> CGA2D {
      let p = CGA2D(e1: q.x, e2: q.y,e: -1,ebar: 0)
        
       let m =   -1 * p * n * p
       
        //X = x^2 n + 2x - ñ,
       let x_sq = length_squared(q)
        print(x_sq)
       let q =  CGA2D(e1: 2 * q.x, e2: 2 * q.y,e: 0,ebar: 0)
        
        let res =  q + x_sq * n - ñ
       
       print(res)
        
        return res
    }
    
    var drop:SIMD2<Double> {
        let v = self.standard
        
        //because the standard form is X = x^2 n + 2x - ñ,
        // we can read of
        let x = 0.5 * v.e1
        let y = 0.5 * v.e2
        
       return [x,y]
    }
    
    var scalerPart:Double {
        values[0]
    }
    
    var vectorPart:Vector {
        Vector(e1: values[1], e2: values[2], e: values[3], ebar: values[4])
    }
    
    var standard:Vector {
        // Standard form of null vector is
        // X = x^2 n + 2x - ñ
        let z = self.vectorPart
        let n = CGA2D.n.vectorPart
        
        let λ = -2 / Vector.dot(z,n)
        
        return (λ * z)
    }
    
    subscript(index: Int) -> Double {
        get {
            values[index]
        }
        set(newValue) {
            values[index] = newValue
        }
    }

}


extension CGA2D {

  //***********************
  // CGA2D.Reverse : res = ~a
  // Reverse the order of the basis blades.
  //***********************

      static prefix func ~ (_ a:CGA2D) -> CGA2D {
      var res = CGA2D.zero
      res[0] = a[0];
  res[1] = a[1];
  res[2] = a[2];
  res[3] = a[3];
  res[4] = a[4];
  res[5] = -a[5];
  res[6] = -a[6];
  res[7] = -a[7];
  res[8] = -a[8];
  res[9] = -a[9];
  res[10] = -a[10];
  res[11] = -a[11];
  res[12] = -a[12];
  res[13] = -a[13];
  res[14] = -a[14];
  res[15] = a[15];
      return res;
    }



  //***********************
  // CGA2D.Dual : res = !a
  // Poincare duality operator.
  //***********************

      static prefix func ! (_ a:CGA2D) -> CGA2D {
      var res = CGA2D.zero
      res[0] = -a[15];
  res[1] = a[14];
  res[2] = -a[13];
  res[3] = a[12];
  res[4] = a[11];
  res[5] = a[10];
  res[6] = -a[9];
  res[7] = -a[8];
  res[8] = a[7];
  res[9] = a[6];
  res[10] = -a[5];
  res[11] = -a[4];
  res[12] = -a[3];
  res[13] = a[2];
  res[14] = -a[1];
  res[15] = a[0];
      return res;
    }



  //***********************
  // CGA2D.Conjugate : res = a.Conjugate()
  // Clifford Conjugation
  //***********************

     func Conjugate (_ a:CGA2D) -> CGA2D {
      var res = CGA2D.zero
      res[0] = a[0];
  res[1] = -a[1];
  res[2] = -a[2];
  res[3] = -a[3];
  res[4] = -a[4];
  res[5] = -a[5];
  res[6] = -a[6];
  res[7] = -a[7];
  res[8] = -a[8];
  res[9] = -a[9];
  res[10] = -a[10];
  res[11] = a[11];
  res[12] = a[12];
  res[13] = a[13];
  res[14] = a[14];
  res[15] = a[15];
      return res;
    }



  //***********************
  // CGA2D.Involute : res = a.Involute()
  // Main involution
  //***********************

     func Involute (_ a:CGA2D) -> CGA2D {
      var res = CGA2D.zero
      res[0] = a[0];
  res[1] = -a[1];
  res[2] = -a[2];
  res[3] = -a[3];
  res[4] = -a[4];
  res[5] = a[5];
  res[6] = a[6];
  res[7] = a[7];
  res[8] = a[8];
  res[9] = a[9];
  res[10] = a[10];
  res[11] = -a[11];
  res[12] = -a[12];
  res[13] = -a[13];
  res[14] = -a[14];
  res[15] = a[15];
      return res;
    }

//***********************
// CGA2D.Mul : res = a * b
// The geometric product.
//***********************

  static func * (_ a:CGA2D, _ b:CGA2D) -> CGA2D{
      var res = CGA2D.zero
      res[0] = b[0] * a[0] + b[1] * a[1] + b[2] * a[2] + b[3] * a[3] - b[4] * a[4] - b[5] * a[5] - b[6] * a[6] + b[7] * a[7] - b[8] * a[8] + b[9] * a[9] + b[10] * a[10] - b[11] * a[11] + b[12] * a[12] + b[13] * a[13] + b[14] * a[14] - b[15] * a[15];
  res[1] = b[1] * a[0] + b[0] * a[1] - b[5] * a[2] - b[6] * a[3] + b[7] * a[4] + b[2] * a[5] + b[3] * a[6] - b[4] * a[7] - b[11] * a[8] + b[12] * a[9] + b[13] * a[10] - b[8] * a[11] + b[9] * a[12] + b[10] * a[13] - b[15] * a[14] + b[14] * a[15];
  res[2] = b[2] * a[0] + b[5] * a[1] + b[0] * a[2] - b[8] * a[3] + b[9] * a[4] - b[1] * a[5] + b[11] * a[6] - b[12] * a[7] + b[3] * a[8] - b[4] * a[9] + b[14] * a[10] + b[6] * a[11] - b[7] * a[12] + b[15] * a[13] + b[10] * a[14] - b[13] * a[15];
  res[3] = b[3] * a[0] + b[6] * a[1] + b[8] * a[2] + b[0] * a[3] + b[10] * a[4] - b[11] * a[5] - b[1] * a[6] - b[13] * a[7] - b[2] * a[8] - b[14] * a[9] - b[4] * a[10] - b[5] * a[11] - b[15] * a[12] - b[7] * a[13] - b[9] * a[14] + b[12] * a[15];
  res[4] = b[4] * a[0] + b[7] * a[1] + b[9] * a[2] + b[10] * a[3] + b[0] * a[4] - b[12] * a[5] - b[13] * a[6] - b[1] * a[7] - b[14] * a[8] - b[2] * a[9] - b[3] * a[10] - b[15] * a[11] - b[5] * a[12] - b[6] * a[13] - b[8] * a[14] + b[11] * a[15];
  res[5] = b[5] * a[0] + b[2] * a[1] - b[1] * a[2] + b[11] * a[3] - b[12] * a[4] + b[0] * a[5] - b[8] * a[6] + b[9] * a[7] + b[6] * a[8] - b[7] * a[9] + b[15] * a[10] + b[3] * a[11] - b[4] * a[12] + b[14] * a[13] - b[13] * a[14] + b[10] * a[15];
  res[6] = b[6] * a[0] + b[3] * a[1] - b[11] * a[2] - b[1] * a[3] - b[13] * a[4] + b[8] * a[5] + b[0] * a[6] + b[10] * a[7] - b[5] * a[8] - b[15] * a[9] - b[7] * a[10] - b[2] * a[11] - b[14] * a[12] - b[4] * a[13] + b[12] * a[14] - b[9] * a[15];
  res[7] = b[7] * a[0] + b[4] * a[1] - b[12] * a[2] - b[13] * a[3] - b[1] * a[4] + b[9] * a[5] + b[10] * a[6] + b[0] * a[7] - b[15] * a[8] - b[5] * a[9] - b[6] * a[10] - b[14] * a[11] - b[2] * a[12] - b[3] * a[13] + b[11] * a[14] - b[8] * a[15];
  res[8] = b[8] * a[0] + b[11] * a[1] + b[3] * a[2] - b[2] * a[3] - b[14] * a[4] - b[6] * a[5] + b[5] * a[6] + b[15] * a[7] + b[0] * a[8] + b[10] * a[9] - b[9] * a[10] + b[1] * a[11] + b[13] * a[12] - b[12] * a[13] - b[4] * a[14] + b[7] * a[15];
  res[9] = b[9] * a[0] + b[12] * a[1] + b[4] * a[2] - b[14] * a[3] - b[2] * a[4] - b[7] * a[5] + b[15] * a[6] + b[5] * a[7] + b[10] * a[8] + b[0] * a[9] - b[8] * a[10] + b[13] * a[11] + b[1] * a[12] - b[11] * a[13] - b[3] * a[14] + b[6] * a[15];
  res[10] = b[10] * a[0] + b[13] * a[1] + b[14] * a[2] + b[4] * a[3] - b[3] * a[4] - b[15] * a[5] - b[7] * a[6] + b[6] * a[7] - b[9] * a[8] + b[8] * a[9] + b[0] * a[10] - b[12] * a[11] + b[11] * a[12] + b[1] * a[13] + b[2] * a[14] - b[5] * a[15];
  res[11] = b[11] * a[0] + b[8] * a[1] - b[6] * a[2] + b[5] * a[3] + b[15] * a[4] + b[3] * a[5] - b[2] * a[6] - b[14] * a[7] + b[1] * a[8] + b[13] * a[9] - b[12] * a[10] + b[0] * a[11] + b[10] * a[12] - b[9] * a[13] + b[7] * a[14] - b[4] * a[15];
  res[12] = b[12] * a[0] + b[9] * a[1] - b[7] * a[2] + b[15] * a[3] + b[5] * a[4] + b[4] * a[5] - b[14] * a[6] - b[2] * a[7] + b[13] * a[8] + b[1] * a[9] - b[11] * a[10] + b[10] * a[11] + b[0] * a[12] - b[8] * a[13] + b[6] * a[14] - b[3] * a[15];
  res[13] = b[13] * a[0] + b[10] * a[1] - b[15] * a[2] - b[7] * a[3] + b[6] * a[4] + b[14] * a[5] + b[4] * a[6] - b[3] * a[7] - b[12] * a[8] + b[11] * a[9] + b[1] * a[10] - b[9] * a[11] + b[8] * a[12] + b[0] * a[13] - b[5] * a[14] + b[2] * a[15];
  res[14] = b[14] * a[0] + b[15] * a[1] + b[10] * a[2] - b[9] * a[3] + b[8] * a[4] - b[13] * a[5] + b[12] * a[6] - b[11] * a[7] + b[4] * a[8] - b[3] * a[9] + b[2] * a[10] + b[7] * a[11] - b[6] * a[12] + b[5] * a[13] + b[0] * a[14] - b[1] * a[15];
  res[15] = b[15] * a[0] + b[14] * a[1] - b[13] * a[2] + b[12] * a[3] - b[11] * a[4] + b[10] * a[5] - b[9] * a[6] + b[8] * a[7] + b[7] * a[8] - b[6] * a[9] + b[5] * a[10] + b[4] * a[11] - b[3] * a[12] + b[2] * a[13] - b[1] * a[14] + b[0] * a[15];
      return res
  }



//***********************
// CGA2D.Wedge : res = a ^ b
// The outer product. (MEET)
//***********************

  static func ^ (_ a:CGA2D, _ b:CGA2D) -> CGA2D{
      var res = CGA2D.zero
      res[0] = b[0] * a[0];
  res[1] = b[1] * a[0] + b[0] * a[1];
  res[2] = b[2] * a[0] + b[0] * a[2];
  res[3] = b[3] * a[0] + b[0] * a[3];
  res[4] = b[4] * a[0] + b[0] * a[4];
  res[5] = b[5] * a[0] + b[2] * a[1] - b[1] * a[2] + b[0] * a[5];
  res[6] = b[6] * a[0] + b[3] * a[1] - b[1] * a[3] + b[0] * a[6];
  res[7] = b[7] * a[0] + b[4] * a[1] - b[1] * a[4] + b[0] * a[7];
  res[8] = b[8] * a[0] + b[3] * a[2] - b[2] * a[3] + b[0] * a[8];
  res[9] = b[9] * a[0] + b[4] * a[2] - b[2] * a[4] + b[0] * a[9];
  res[10] = b[10] * a[0] + b[4] * a[3] - b[3] * a[4] + b[0] * a[10];
  res[11] = b[11] * a[0] + b[8] * a[1] - b[6] * a[2] + b[5] * a[3] + b[3] * a[5] - b[2] * a[6] + b[1] * a[8] + b[0] * a[11];
  res[12] = b[12] * a[0] + b[9] * a[1] - b[7] * a[2] + b[5] * a[4] + b[4] * a[5] - b[2] * a[7] + b[1] * a[9] + b[0] * a[12];
  res[13] = b[13] * a[0] + b[10] * a[1] - b[7] * a[3] + b[6] * a[4] + b[4] * a[6] - b[3] * a[7] + b[1] * a[10] + b[0] * a[13];
  res[14] = b[14] * a[0] + b[10] * a[2] - b[9] * a[3] + b[8] * a[4] + b[4] * a[8] - b[3] * a[9] + b[2] * a[10] + b[0] * a[14];
  res[15] = b[15] * a[0] + b[14] * a[1] - b[13] * a[2] + b[12] * a[3] - b[11] * a[4] + b[10] * a[5] - b[9] * a[6] + b[8] * a[7] + b[7] * a[8] - b[6] * a[9] + b[5] * a[10] + b[4] * a[11] - b[3] * a[12] + b[2] * a[13] - b[1] * a[14] + b[0] * a[15];
      return res
  }



//***********************
// CGA2D.Vee : res = a & b
// The regressive product. (JOIN)
//***********************

  static func & (_ a:CGA2D, _ b:CGA2D) -> CGA2D{
      var res = CGA2D.zero
      res[15] = 1 * (a[15] * b[15]);
      res[14] =  -1 * (a[14] *  -1 * b[15] + a[15] * b[14] *  -1);
  res[13] = 1 * (a[13] * b[15] + a[15] * b[13]);
      res[12] =  -1 * (a[12] *  -1 * b[15] + a[15] * b[12] *  -1);
  res[11] = 1 * (a[11] * b[15] + a[15] * b[11]);
      res[10] = 1 * (a[10] * b[15] + a[13] * b[14] *  -1 - a[14] *  -1 * b[13] + a[15] * b[10]);
      res[9] =  -1 * (a[9] *  -1 * b[15] + a[12] *  -1 * b[14] *  -1 - a[14] *  -1 * b[12] *  -1 + a[15] * b[9] *  -1);
      res[8] = 1 * (a[8] * b[15] + a[11] * b[14] *  -1 - a[14] *  -1 * b[11] + a[15] * b[8]);
      res[7] = 1 * (a[7] * b[15] + a[12] *  -1 * b[13] - a[13] * b[12] *  -1 + a[15] * b[7]);
      res[6] =  -1 * (a[6] *  -1 * b[15] + a[11] * b[13] - a[13] * b[11] + a[15] * b[6] *  -1);
      res[5] = 1 * (a[5] * b[15] + a[11] * b[12] *  -1 - a[12] *  -1 * b[11] + a[15] * b[5]);
      res[4] =  -1 * (a[4] *  -1 * b[15] + a[7] * b[14] *  -1 - a[9] *  -1 * b[13] + a[10] * b[12] *  -1 + a[12] *  -1 * b[10] - a[13] * b[9] *  -1 + a[14] *  -1 * b[7] + a[15] * b[4] *  -1);
      res[3] = 1 * (a[3] * b[15] + a[6] *  -1 * b[14] *  -1 - a[8] * b[13] + a[10] * b[11] + a[11] * b[10] - a[13] * b[8] + a[14] *  -1 * b[6] *  -1 + a[15] * b[3]);
      res[2] =  -1 * (a[2] *  -1 * b[15] + a[5] * b[14] *  -1 - a[8] * b[12] *  -1 + a[9] *  -1 * b[11] + a[11] * b[9] *  -1 - a[12] *  -1 * b[8] + a[14] *  -1 * b[5] + a[15] * b[2] *  -1);
      res[1] = 1 * (a[1] * b[15] + a[5] * b[13] - a[6] *  -1 * b[12] *  -1 + a[7] * b[11] + a[11] * b[7] - a[12] *  -1 * b[6] *  -1 + a[13] * b[5] + a[15] * b[1]);
      res[0] = 1 * (a[0] * b[15] + a[1] * b[14] *  -1 - a[2] *  -1 * b[13] + a[3] * b[12] *  -1 - a[4] *  -1 * b[11] + a[5] * b[10] - a[6] *  -1 * b[9] *  -1 + a[7] * b[8] + a[8] * b[7] - a[9] *  -1 * b[6] *  -1 + a[10] * b[5] + a[11] * b[4] *  -1 - a[12] *  -1 * b[3] + a[13] * b[2] *  -1 - a[14] *  -1 * b[1] + a[15] * b[0]);
      return res
  }



//***********************
// CGA2D.Dot : res = a | b
// The inner product.
//***********************

  static func | (_ a:CGA2D, _ b:CGA2D) -> CGA2D {
      var res = CGA2D.zero
      res[0] = b[0] * a[0] + b[1] * a[1] + b[2] * a[2] + b[3] * a[3] - b[4] * a[4] - b[5] * a[5] - b[6] * a[6] + b[7] * a[7] - b[8] * a[8] + b[9] * a[9] + b[10] * a[10] - b[11] * a[11] + b[12] * a[12] + b[13] * a[13] + b[14] * a[14] - b[15] * a[15];
  res[1] = b[1] * a[0] + b[0] * a[1] - b[5] * a[2] - b[6] * a[3] + b[7] * a[4] + b[2] * a[5] + b[3] * a[6] - b[4] * a[7] - b[11] * a[8] + b[12] * a[9] + b[13] * a[10] - b[8] * a[11] + b[9] * a[12] + b[10] * a[13] - b[15] * a[14] + b[14] * a[15];
  res[2] = b[2] * a[0] + b[5] * a[1] + b[0] * a[2] - b[8] * a[3] + b[9] * a[4] - b[1] * a[5] + b[11] * a[6] - b[12] * a[7] + b[3] * a[8] - b[4] * a[9] + b[14] * a[10] + b[6] * a[11] - b[7] * a[12] + b[15] * a[13] + b[10] * a[14] - b[13] * a[15];
  res[3] = b[3] * a[0] + b[6] * a[1] + b[8] * a[2] + b[0] * a[3] + b[10] * a[4] - b[11] * a[5] - b[1] * a[6] - b[13] * a[7] - b[2] * a[8] - b[14] * a[9] - b[4] * a[10] - b[5] * a[11] - b[15] * a[12] - b[7] * a[13] - b[9] * a[14] + b[12] * a[15];
  res[4] = b[4] * a[0] + b[7] * a[1] + b[9] * a[2] + b[10] * a[3] + b[0] * a[4] - b[12] * a[5] - b[13] * a[6] - b[1] * a[7] - b[14] * a[8] - b[2] * a[9] - b[3] * a[10] - b[15] * a[11] - b[5] * a[12] - b[6] * a[13] - b[8] * a[14] + b[11] * a[15];
  res[5] = b[5] * a[0] + b[11] * a[3] - b[12] * a[4] + b[0] * a[5] + b[15] * a[10] + b[3] * a[11] - b[4] * a[12] + b[10] * a[15];
  res[6] = b[6] * a[0] - b[11] * a[2] - b[13] * a[4] + b[0] * a[6] - b[15] * a[9] - b[2] * a[11] - b[4] * a[13] - b[9] * a[15];
  res[7] = b[7] * a[0] - b[12] * a[2] - b[13] * a[3] + b[0] * a[7] - b[15] * a[8] - b[2] * a[12] - b[3] * a[13] - b[8] * a[15];
  res[8] = b[8] * a[0] + b[11] * a[1] - b[14] * a[4] + b[15] * a[7] + b[0] * a[8] + b[1] * a[11] - b[4] * a[14] + b[7] * a[15];
  res[9] = b[9] * a[0] + b[12] * a[1] - b[14] * a[3] + b[15] * a[6] + b[0] * a[9] + b[1] * a[12] - b[3] * a[14] + b[6] * a[15];
  res[10] = b[10] * a[0] + b[13] * a[1] + b[14] * a[2] - b[15] * a[5] + b[0] * a[10] + b[1] * a[13] + b[2] * a[14] - b[5] * a[15];
  res[11] = b[11] * a[0] + b[15] * a[4] + b[0] * a[11] - b[4] * a[15];
  res[12] = b[12] * a[0] + b[15] * a[3] + b[0] * a[12] - b[3] * a[15];
  res[13] = b[13] * a[0] - b[15] * a[2] + b[0] * a[13] + b[2] * a[15];
  res[14] = b[14] * a[0] + b[15] * a[1] + b[0] * a[14] - b[1] * a[15];
  res[15] = b[15] * a[0] + b[0] * a[15];
      return res
  }



//***********************
// CGA2D.Add : res = a + b
// Multivector addition
//***********************

  static func + (_ a:CGA2D, _ b:CGA2D) -> CGA2D{
      var res = CGA2D.zero
          res[0]  =  a[0] + b[0];
    res[1]  =  a[1] + b[1];
    res[2]  =  a[2] + b[2];
    res[3]  =  a[3] + b[3];
    res[4]  =  a[4] + b[4];
    res[5]  =  a[5] + b[5];
    res[6]  =  a[6] + b[6];
    res[7]  =  a[7] + b[7];
    res[8]  =  a[8] + b[8];
    res[9]  =  a[9] + b[9];
    res[10]  =  a[10] + b[10];
    res[11]  =  a[11] + b[11];
    res[12]  =  a[12] + b[12];
    res[13]  =  a[13] + b[13];
    res[14]  =  a[14] + b[14];
    res[15]  =  a[15] + b[15];
      return res
  }



//***********************
// CGA2D.Sub : res = a - b
// Multivector subtraction
//***********************

  static func - (_ a:CGA2D, _ b:CGA2D) -> CGA2D{
      var res = CGA2D.zero
          res[0]  =  a[0] - b[0];
    res[1]  =  a[1] - b[1];
    res[2]  =  a[2] - b[2];
    res[3]  =  a[3] - b[3];
    res[4]  =  a[4] - b[4];
    res[5]  =  a[5] - b[5];
    res[6]  =  a[6] - b[6];
    res[7]  =  a[7] - b[7];
    res[8]  =  a[8] - b[8];
    res[9]  =  a[9] - b[9];
    res[10]  =  a[10] - b[10];
    res[11]  =  a[11] - b[11];
    res[12]  =  a[12] - b[12];
    res[13]  =  a[13] - b[13];
    res[14]  =  a[14] - b[14];
    res[15]  =  a[15] - b[15];
      return res
  }



//***********************
// CGA2D.smul : res = a * b
// scalar/multivector multiplication
//***********************

  static func * (_ a:Double, _ b:CGA2D) -> CGA2D{
      var res = CGA2D.zero
          res[0]  =  a * b[0];
    res[1]  =  a * b[1];
    res[2]  =  a * b[2];
    res[3]  =  a * b[3];
    res[4]  =  a * b[4];
    res[5]  =  a * b[5];
    res[6]  =  a * b[6];
    res[7]  =  a * b[7];
    res[8]  =  a * b[8];
    res[9]  =  a * b[9];
    res[10]  =  a * b[10];
    res[11]  =  a * b[11];
    res[12]  =  a * b[12];
    res[13]  =  a * b[13];
    res[14]  =  a * b[14];
    res[15]  =  a * b[15];
      return res
  }



//***********************
// CGA2D.muls : res = a * b
// multivector/scalar multiplication
//***********************

  static func * (_ a:CGA2D, _ b:Double) -> CGA2D{
      var res = CGA2D.zero
          res[0]  =  a[0] * b;
    res[1]  =  a[1] * b;
    res[2]  =  a[2] * b;
    res[3]  =  a[3] * b;
    res[4]  =  a[4] * b;
    res[5]  =  a[5] * b;
    res[6]  =  a[6] * b;
    res[7]  =  a[7] * b;
    res[8]  =  a[8] * b;
    res[9]  =  a[9] * b;
    res[10]  =  a[10] * b;
    res[11]  =  a[11] * b;
    res[12]  =  a[12] * b;
    res[13]  =  a[13] * b;
    res[14]  =  a[14] * b;
    res[15]  =  a[15] * b;
      return res
  }



//***********************
// CGA2D.sadd : res = a + b
// scalar/multivector addition
//***********************

  static func + (_ a:Double, _ b:CGA2D) -> CGA2D{
      var res = CGA2D.zero
        res[0]  =  a + b[0];
      res[1]  =  b[1];
    res[2]  =  b[2];
    res[3]  =  b[3];
    res[4]  =  b[4];
    res[5]  =  b[5];
    res[6]  =  b[6];
    res[7]  =  b[7];
    res[8]  =  b[8];
    res[9]  =  b[9];
    res[10]  =  b[10];
    res[11]  =  b[11];
    res[12]  =  b[12];
    res[13]  =  b[13];
    res[14]  =  b[14];
    res[15]  =  b[15];
      return res
  }



//***********************
// CGA2D.adds : res = a + b
// multivector/scalar addition
//***********************

  static func + (_ a:CGA2D, _ b:Double) -> CGA2D{
      var res = CGA2D.zero
        res[0]  =  a[0] + b;
      res[1]  =  a[1];
    res[2]  =  a[2];
    res[3]  =  a[3];
    res[4]  =  a[4];
    res[5]  =  a[5];
    res[6]  =  a[6];
    res[7]  =  a[7];
    res[8]  =  a[8];
    res[9]  =  a[9];
    res[10]  =  a[10];
    res[11]  =  a[11];
    res[12]  =  a[12];
    res[13]  =  a[13];
    res[14]  =  a[14];
    res[15]  =  a[15];
      return res
  }



//***********************
// CGA2D.ssub : res = a - b
// scalar/multivector subtraction
//***********************

  static func - (_ a:Double, _ b:CGA2D) -> CGA2D{
      var res = CGA2D.zero
        res[0]  =  a - b[0];
      res[1]  =   -b[1];
      res[2]  =   -b[2];
      res[3]  =   -b[3];
      res[4]  =   -b[4];
      res[5]  =   -b[5];
      res[6]  =   -b[6];
      res[7]  =   -b[7];
      res[8]  =   -b[8];
      res[9]  =   -b[9];
      res[10]  =   -b[10];
      res[11]  =   -b[11];
      res[12]  =   -b[12];
      res[13]  =   -b[13];
      res[14]  =   -b[14];
      res[15]  =   -b[15];
      return res
  }



//***********************
// CGA2D.subs : res = a - b
// multivector/scalar subtraction
//***********************

  static func - (_ a:CGA2D, _ b:Double) -> CGA2D{
      var res = CGA2D.zero
        res[0]  =  a[0] - b;
      res[1]  =  a[1];
    res[2]  =  a[2];
    res[3]  =  a[3];
    res[4]  =  a[4];
    res[5]  =  a[5];
    res[6]  =  a[6];
    res[7]  =  a[7];
    res[8]  =  a[8];
    res[9]  =  a[9];
    res[10]  =  a[10];
    res[11]  =  a[11];
    res[12]  =  a[12];
    res[13]  =  a[13];
    res[14]  =  a[14];
    res[15]  =  a[15];
      return res
  }
   
    public static func circle(x:CGA2D,y:CGA2D,z:CGA2D) -> CGA2D {
        x ^ y ^ z
    }
    
    public func circleCenterAndRadius() -> (center:SIMD2<Double>,radius:Double) {
        let L = self
        
        let Lsq = (L | L).scalerPart
        
        let  L_wedge_n = (L ^ CGA2D.n)
        
        let L_wedge_n_sq = (L_wedge_n * L_wedge_n).scalerPart
        
        print( "L^2 \(Lsq). L_wedge_n^2 scaler part: \(L_wedge_n_sq)")
        
        let r_squard = -Lsq / L_wedge_n_sq
        // see Geometric Algebra for Physists eqn 10.106
        
        let C = L * CGA2D.n * L
        
        return (center:C.drop,radius:sqrt(r_squard))
    }
    
    var isNull:Bool { norm.isApproximatelyEqual(to: 0,absoluteTolerance: 0.0001) }
    
    var norm:Double {

        let  res = self * Conjugate(self)
        let scalar = res[0]
        
        return sqrt(abs(scalar))
    }

    var inorm:Double {
        let dual = !(self);
        return dual.norm;
    }

    var normalized:CGA2D {
        self * (1/norm);
    }

    func log() {
        let a = self
        print("Scaler \(a[0])")
        print("Vector \(a[1]),\(a[2]),\(a[3]),\(a[4])")
        print("BiVector \(a[5]),\(a[6]),\(a[7]),\(a[8]),\(a[9]),\(a[10])")
        print("TriVector \(a[11]),\(a[12]),\(a[13]),\(a[14]))")
        print("PuedoScaler \(a[15])")
    }
}

extension CGA2D : CustomStringConvertible {
    public var description: String {
        let a = self
        return  """
                Scaler \(a[0])
                Vector \(a[1]),\(a[2]),\(a[3]),\(a[4])
                BiVector \(a[5]),\(a[6]),\(a[7]),\(a[8]),\(a[9]),\(a[10])
                TriVector \(a[11]),\(a[12]),\(a[13]),\(a[14]))
                PuedoScaler \(a[15])
                """
    }
    
    
}



fileprivate func example() {
  
  //print("e1*e1         : "); (e1*e1).log();
  //print("pss           : "); e1234.log();
  //print("pss*pss       : "); (e1234*e1234).log();

}



