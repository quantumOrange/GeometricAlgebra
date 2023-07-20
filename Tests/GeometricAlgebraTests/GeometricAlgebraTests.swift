import XCTest
@testable import GeometricAlgebra
let accuracy:Double = 0.0001

final class GeometricAlgebraTests: XCTestCase {
    
    
    func testConvertFromPlane() throws {
        let o = SIMD2<Double>(0,0)
        let p = SIMD2<Double>(0.21,0.21)
        let q = SIMD2<Double>(-2.0,-4)
        let s = SIMD2<Double>(0.4,1.7)
        
        let O =  34 * CGA2D.lift(point: o)
        let P =  34 * CGA2D.lift(point: p)
        let Q = 0.3 * CGA2D.lift(point: q)
        let S = 0.7 * CGA2D.lift(point: s)
        
        
       XCTAssertEqual(O.norm, 0, accuracy: accuracy)
       XCTAssertEqual(P.norm, 0, accuracy: accuracy)
       XCTAssertEqual(Q.norm, 0, accuracy: accuracy)
       XCTAssertEqual(S.norm, 0, accuracy: accuracy)
        
        XCTAssert(P.isNull)
        XCTAssert(Q.isNull)
        XCTAssert(S.isNull)
        
        
        let p_ = P.drop
        let q_ = Q.drop
        let s_ = S.drop
        

        XCTAssertEqual(p.x,p_.x,accuracy: accuracy)
        XCTAssertEqual(p.y,p_.y,accuracy: accuracy)
        XCTAssertEqual(q.x,q_.x,accuracy: accuracy)
        XCTAssertEqual(q.y,q_.y,accuracy: accuracy)
        XCTAssertEqual(s.x,s_.x,accuracy: accuracy)
        XCTAssertEqual(s.y,s_.y,accuracy: accuracy)
 
    }
    
    func testConvertFromPlane2() throws {
       
        let p = SIMD2<Double>(0.5,0.03)
        let q = SIMD2<Double>(-0.1,-3)
        let s = SIMD2<Double>(0.2,-0.12)
        
        
        let P =  34 * CGA2D.lift(point: p)
        let Q = 0.3 * CGA2D.lift(point: q)
        let S = 0.7 * CGA2D.lift(point: s)
        
        XCTAssert(P.isNull)
        XCTAssert(Q.isNull)
        XCTAssert(S.isNull)
        
        let p_ = P.drop
        let q_ = Q.drop
        let s_ = S.drop
        
        XCTAssertEqual(p.x,p_.x,accuracy: accuracy)
        XCTAssertEqual(p.y,p_.y,accuracy: accuracy)
        XCTAssertEqual(q.x,q_.x,accuracy: accuracy)
        XCTAssertEqual(q.y,q_.y,accuracy: accuracy)
        XCTAssertEqual(s.x,s_.x,accuracy: accuracy)
        XCTAssertEqual(s.y,s_.y,accuracy: accuracy)
 
    }
    
    func testCircle() throws {
        // These three points define a circle with radius 1 and center ant the origin (0,0)
        let p = SIMD2<Double>(1,0)
        let q = SIMD2<Double>(-1,0)
        let s = SIMD2<Double>(0,1)
        
        let P = CGA2D.lift(point: p)
        let Q = CGA2D.lift(point: q)
        let S = CGA2D.lift(point: s)
        
        let L = P ^ Q ^ S
        
        let (center,radius) = L.circleCenterAndRadius()
        
        XCTAssertEqual(center.x, 0)
        XCTAssertEqual(center.y, 0)
        XCTAssertEqual(radius, 1)
        
    }
    
    func testSmallCircle() throws {
        // These three points define a circle with radius 0.5 and center ant the origin (0,0)
        let p = SIMD2<Double>(0.5,0)
        let q = SIMD2<Double>(-0.5,0)
        let s = SIMD2<Double>(0,0.5)
        
        let P = CGA2D.lift(point: p)
        let Q = CGA2D.lift(point: q)
        let S = CGA2D.lift(point: s)
        
        let L = P ^ Q ^ S
        
        let (center,radius) = L.circleCenterAndRadius()
        
        XCTAssertEqual(center.x, 0)
        XCTAssertEqual(center.y, 0)
        XCTAssertEqual(radius, 0.5)
        
    }
    
    func testCircle3() throws {
        // These three points define a circle with radius 0.5 and center ant the origin (1,1)
        let t = SIMD2<Double>(1,1)
        let p = SIMD2<Double>(0.5,0) + t
        let q = SIMD2<Double>(-0.5,0)  + t
        let s = SIMD2<Double>(0,0.5)  + t
        
        
         
        let P = CGA2D.lift(point: p)
        let Q = CGA2D.lift(point: q)
        let S = CGA2D.lift(point: s)
        
        let L = P ^ Q ^ S
        
        let (center,radius) = L.circleCenterAndRadius()
        
        XCTAssertEqual(center.x, 1)
        XCTAssertEqual(center.y, 1)
        XCTAssertEqual(radius, 0.5)
        
    }
    
    
    /*
    func testDistance() {
        
        let x = SIMD2<Double>(-0.1,0.3)
        let y = SIMD2<Double>(0.23,0.3)
        
        let d1 = hyperbolicDistance2(p: x,q: y)
        
        let X = HyperCga2D.lift(point: x)
        let Y = HyperCga2D.lift(point: y)
        
        let d2 = HyperCga2D.hyperbolicDistance(x: X,y: Y)
        
        XCTAssertEqual(d1,d2,accuracy: accuracy)
        
    }
     
    
    func testTranslation() {
        
        let x = SIMD2<Double>(0.12,0.53)
        let y = SIMD2<Double>(0.5,0.1)
        
        let a = SIMD2<Double>(0.1,0.1)
        let b = SIMD2<Double>(0.05,0.25)
        
        
        
    }
     
     */
}
