//
//  File.swift
//  
//
//  Created by David Crooks on 15/07/2023.
//

import Foundation
import simd
import Numerics

extension HyperCga2D {
   
    public static func lift(point q:SIMD2<Double>) -> HyperCga2D {
      let p = HyperCga2D(e1: q.x, e2: q.y,e: -1,ebar: 0)
        
    let m =   -1 * p * n * p
        
        //X = x^2 n + 2x - ñ,
       let x_sq = length_squared(q)
        print(x_sq)
       let q =  HyperCga2D(e1: 2 * q.x, e2: 2 * q.y,e: 0,ebar: 0)
        
        let res =  q + x_sq * n - ñ
       
       print(res)
        
        return res
    }
    
    public var drop:SIMD2<Double> {
        let v = self.standard
        
        //because the standard form is X = x^2 n + 2x - ñ,
        // we can read of
        let x = 0.5 * v.e1
        let y = 0.5 * v.e2
        
       return [x,y]
    }
    
   
    
    var standard:Vector {
        // Standard form of null vector is
        // X = x^2 n + 2x - ñ
        let z = self.vectorPart
        let n = HyperCga2D.n.vectorPart
        
        let λ = -2 / Vector.dot(z,n)
        
        return (λ * z)
    }
    
   
    public static func dline(x:HyperCga2D,y:HyperCga2D) -> HyperCga2D {
        x ^ y ^ e
    }
    
    public static func translation(x:HyperCga2D,y:HyperCga2D) -> HyperCga2D {
        (x ^ y ^ e) * e
    }
    
    public static func scaledTranslation(x:HyperCga2D,y:HyperCga2D) -> HyperCga2D {
        let T = translation(x: x, y: y)
        let d1 = hyperbolicDistance(x: x, y: y)
        let y_dash = x.apply(T)
        let d2 = hyperbolicDistance(x: x, y: y_dash)
        guard d2 > 0.0000001 else  { return .zero }
        return (d1 / d2) * T
    }
    
    
    public static func exp(_ a:Double = 1.0, _ B:HyperCga2D) -> HyperCga2D {
        let k = a * B.norm
        let B_hat = B.normalized
        let B_sq = (B * B).scalerPart
        
        if B_sq > 0 {
            return cosh(k) + sinh(k) * B_hat
        }
        else if B_sq < 0 {
            return cos(k) + sin(k) * B_hat
        }
        return 1 + k * B_hat
    }
    
    
    public static  func shouldBeScalar(_ T:HyperCga2D) -> HyperCga2D{
        T * ~T
    }
    
    public var sq:HyperCga2D{
        self * self
    }
    
    public  func applyRotor(_ V:HyperCga2D) -> HyperCga2D{
        V * self * ~V
    }
    
    public  func apply(_ V:HyperCga2D, alpha α:Double = 1.0) -> HyperCga2D {
        HyperCga2D.exp( α/2, V) * self * HyperCga2D.exp( -α/2, V)
    }
    
    
    public static func hyperbolicDistance(x:HyperCga2D,y:HyperCga2D)  -> Double {
        hyperbolicDistance(x:x.vectorPart,y:y.vectorPart)
    }
    
    static func hyperbolicDistance(x:Vector,y:Vector)  -> Double{
        let numerator = -Vector.dot(x,y)
        
        let x_dot_e = Vector.dot(x,e.vectorPart)
        let y_dot_e = Vector.dot(y,e.vectorPart)

        let denominator =  2 * x_dot_e * y_dot_e
        
        let theta = sqrt(numerator / denominator)
        //print("**------------**")
        //print("cga dist:")
        //print("theta: \(theta) asinsh: \(asinh(theta))  num:\(numerator)  denom: \(denominator) ")
        return 2 * asinh(theta)
    }
}

public func hyperbolicDistance(x:SIMD2<Double>,y:SIMD2<Double>) -> Double {
    let x_sq = dot(x,x)
    let y_sq = dot(y,y)
    
    let numerator = dot(x - y,x - y)
    let denominator = (( 1 - x_sq ) * ( 1 - y_sq ))
    let theta = sqrt(numerator / denominator)
    //print("**------------**")
    //print("simd dist:")
    //print("theta: \(theta) asinsh: \(asinh(theta))  num:\(numerator)  denom: \(denominator) ")
    return 2 * asinh(theta)
}

public  func hyperbolicDistance2(p:SIMD2<Double>, q:SIMD2<Double>) -> Double {
    //distance between points on  Poincaré Hyperbolic Disk
    let pq = length(p-q)
    let op = length(p)
    let oq = length(q)
    
    return acosh(1.0 + 2.0*pq*pq/((1.0 - op*op)*(1.0 - oq*oq)) )
}


func liftToHyperbolid(point q:SIMD2<Double>) -> SIMD3<Double> {
    let r = length(q)
    
 //   let R = -2 * r / (r * r - 1)
    let R = -2  / (r * r - 1)
    let t = (1 +  r * r) /  (1 - r * r)
   
    
    let p = SIMD3<Double>( R * q.x ,R * q.y, t)
    print("hypebolid \(p.x * p.x + p.y * p.y - p.z * p.z)")
    
    
    return p
}

func dropFromHyperbolid(point q:SIMD3<Double>) -> SIMD2<Double> {
   
    let t = q.z
    
    let r_sq  = (t - 1) / (t + 1)
    let R_inv = -(r_sq - 1) / 2
    
    let x = R_inv * q.x
    let y = R_inv * q.y
    
    return [x,y]
}


extension HyperCga2D {
    public var matrix:matrix_float4x4 {
        var m = matrix_float4x4()

        for i in 0..<16 {
            m[i / 4][i % 4] = Float(values[i])
        }
        
        return m
    }
}
