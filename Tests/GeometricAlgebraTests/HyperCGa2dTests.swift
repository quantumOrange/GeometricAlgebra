//
//  HyperCGa2dTests.swift
//  
//
//  Created by David Crooks on 19/07/2023.
//

import XCTest
@testable import GeometricAlgebra
import simd

final class HyperCGa2dTests: XCTestCase {

    override func setUpWithError() throws {
        // Put setup code here. This method is called before the invocation of each test method in the class.
    }

    override func tearDownWithError() throws {
        // Put teardown code here. This method is called after the invocation of each test method in the class.
    }

    
    
    
    func testConvertFromPlane() throws {

        let p = SIMD2<Double>(0.21,0.21)
        let q = SIMD2<Double>(-0.3,-0.21)
        let s = SIMD2<Double>(0.4,0.13)
        // Tests are only valid if these points are in the disk
        assert(length(p) < 1.0)
        assert(length(q) < 1.0)
        assert(length(s) < 1.0)
        
        let P = HyperCga2D.lift(point: p)
        let Q = HyperCga2D.lift(point: q)
        let S = HyperCga2D.lift(point: s)
        
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
    
    func testConvertFromPlaneWithScale() throws {
       
        let p = SIMD2<Double>(0.5,0.03)
        let q = SIMD2<Double>(-0.1,-0.5)
        let s = SIMD2<Double>(0.2,-0.12)
        
        // Tests are only valid if these points are in the disk
        assert(length(p) < 1.0)
        assert(length(q) < 1.0)
        assert(length(s) < 1.0)
        
        let P =  34 * HyperCga2D.lift(point: p)
        let Q = 0.3 * HyperCga2D.lift(point: q)
        let S = 0.7 * HyperCga2D.lift(point: s)
        let S_2 = -0.3 * HyperCga2D.lift(point: s)
        
        XCTAssert(P.isNull)
        XCTAssert(Q.isNull)
        XCTAssert(S.isNull)
        
        let p_ = P.drop
        let q_ = Q.drop
        let s_ = S.drop
        let s_2 = S_2.drop
        
        XCTAssertEqual(p.x,p_.x,accuracy: accuracy)
        XCTAssertEqual(p.y,p_.y,accuracy: accuracy)
        XCTAssertEqual(q.x,q_.x,accuracy: accuracy)
        XCTAssertEqual(q.y,q_.y,accuracy: accuracy)
        XCTAssertEqual(s.x,s_.x,accuracy: accuracy)
        XCTAssertEqual(s.y,s_.y,accuracy: accuracy)
        XCTAssertEqual(s.x,s_2.x,accuracy: accuracy)
        XCTAssertEqual(s.y,s_2.y,accuracy: accuracy)
 
    }
    /*
    func testCircle() throws {
        // These three points define a circle with radius 1 and center ant the origin (0,0)
        let p = SIMD2<Double>(1,0)
        let q = SIMD2<Double>(-1,0)
        let s = SIMD2<Double>(0,1)
        
        let P = HyperCga2D.lift(point: p)
        let Q = HyperCga2D.lift(point: q)
        let S = HyperCga2D.lift(point: s)
        
        let L = P ^ Q ^ S
        
        let (center,radius) = L.circleCenterAndRadius()
        
        XCTAssertEqual(center.x, 0)
        XCTAssertEqual(center.y, 0)
        XCTAssertEqual(radius, 1)
        
    }
    */
    func testSmallCircle() throws {
        // These three points define a circle with radius 0.5 and center ant the origin (0,0)
        let p = SIMD2<Double>(0.5,0)
        let q = SIMD2<Double>(-0.5,0)
        let s = SIMD2<Double>(0,0.5)
        let hyp = hyperbolicDistance2(p: .zero, q: p)
        let P = HyperCga2D.lift(point: p)
        let Q = HyperCga2D.lift(point: q)
        let S = HyperCga2D.lift(point: s)
        
        let L = P ^ Q ^ S
        
        let (center,radius) = L.circleCenterAndRadius()
        let (centerE,radiusE) = L.circleCenterAndRadiusE()
        
        XCTAssertEqual(center.x, 0)
        XCTAssertEqual(center.y, 0)
        XCTAssertEqual(radius, 0.5)
        
        XCTAssertEqual(centerE.x, 0)
        XCTAssertEqual(centerE.y, 0)
        XCTAssertEqual(radiusE, hyp)
        
    }
    
    func testCircle3() throws {
        // These three points define a circle with radius 0.5 and center ant the origin (1,1)
        let t = SIMD2<Double>(1,1)
        let p = SIMD2<Double>(0.5,0) + t
        let q = SIMD2<Double>(-0.5,0)  + t
        let s = SIMD2<Double>(0,0.5)  + t
        
        
         
        let P = HyperCga2D.lift(point: p)
        let Q = HyperCga2D.lift(point: q)
        let S = HyperCga2D.lift(point: s)
        
        let L = P ^ Q ^ S
        
        let (center,radius) = L.circleCenterAndRadius()
        
        XCTAssertEqual(center.x, 1)
        XCTAssertEqual(center.y, 1)
        XCTAssertEqual(radius, 0.5)
        
    }
    
    
    
    func testDistance() {
        
        let p = SIMD2<Double>(-0.1,0.3)
        let q = SIMD2<Double>(0.23,0.3)
        
        let d1 = hyperbolicDistance2(p: p,q: q)
        
        let P = HyperCga2D.lift(point: p)
        let Q = HyperCga2D.lift(point: q)
        
        let d2 = HyperCga2D.hyperbolicDistance(x: P,y: Q)
        
        XCTAssertEqual(d1,d2,accuracy: accuracy)
        
    }
    
    func testExp2() {
        let x = SIMD2<Double>(-0.1,0.3)
        let y = SIMD2<Double>(0.23,0.3)
        let X = HyperCga2D.lift(point: x)
        let Y = HyperCga2D.lift(point: y)
        let e = HyperCga2D.e
        let B =   ( X ^ Y ^ e) * e
        let a = 0.345
        let C = a * B
        let k = B.norm
        let j = C.norm
        
        XCTAssertEqual(a * k,j,accuracy: accuracy)
    }
    
    func testExp() {
        let x = SIMD2<Double>(-0.1,0.3)
        let y = SIMD2<Double>(0.23,0.3)
        let X = HyperCga2D.lift(point: x)
        let Y = HyperCga2D.lift(point: y)
        let e = HyperCga2D.e
        let B =   ( X ^ Y ^ e) * e
        
        let b_sq = (B * B)
        XCTAssert(b_sq.isScaler)
        XCTAssert(b_sq.scalerPart > 0)
        
        let k = B.norm
        let B_hat = B.normalized
        
        XCTAssertEqual((B_hat * B_hat).scalerPart ,1,accuracy: 0.0001)
        
        
        let B_dash = k * B_hat
        //print(B_dash,)
        print(B_dash)
        print("----")
        print(B)
        XCTAssert(B_dash.isApproxmatlyEqual(to: B))
        
    }
    
    
    func testTranslation() {
        
        let p = SIMD2<Double>(-0.1,0.3)
        let q = SIMD2<Double>(0.23,0.3)
        
        let a = SIMD2<Double>(0.1,0.1)
        let b = SIMD2<Double>(0.05,0.25)
        
        let P = HyperCga2D.lift(point: p)
        let Q = HyperCga2D.lift(point: q)
        let A = HyperCga2D.lift(point: a)
        let B = HyperCga2D.lift(point: b)
        
        let T = HyperCga2D.translation(x: P, y: Q)
        
        let P_dash = P.apply(T)
        let A_dash = A.apply(T)
        let B_dash = B.apply(T)
        
        //
        XCTAssertEqual(HyperCga2D.hyperbolicDistance(x: P_dash, y: Q) , 0.0,accuracy: accuracy)
        let p_dash = P_dash.drop
        XCTAssertEqual(p_dash.x,q.x,accuracy: accuracy)
        XCTAssertEqual(p_dash.y,q.y,accuracy: accuracy)
        
        // If T is a translation, then we expect the distances PA, PB and AB to be invarient:
       
        XCTAssertEqual(HyperCga2D.hyperbolicDistance(x: P_dash, y: A_dash),
                       HyperCga2D.hyperbolicDistance(x: P, y: A),accuracy: accuracy)
        
        XCTAssertEqual(HyperCga2D.hyperbolicDistance(x: P_dash, y: B_dash),
                       HyperCga2D.hyperbolicDistance(x: P, y: B),accuracy: accuracy)
        
        XCTAssertEqual(HyperCga2D.hyperbolicDistance(x: A_dash, y: B_dash),
                       HyperCga2D.hyperbolicDistance(x: A, y: B),accuracy: accuracy)
        
        
        // And back again  - testing that ~T inverts the translation
        let P_2 = P_dash.apply(~T)
        let A_2 = A_dash.apply(~T)
        let B_2 = B_dash.apply(~T)
        
        XCTAssertEqual(HyperCga2D.hyperbolicDistance(x: P, y: P_2),
                                                     0.0,accuracy: accuracy)
        
        
        //Thes next 2 tests will fail if the
        XCTAssertEqual(HyperCga2D.hyperbolicDistance(x: A.vectorPart, y: A_2.standard),
                                                     0.0,accuracy: accuracy)
        
        XCTAssertEqual(HyperCga2D.hyperbolicDistance(x: B.vectorPart, y: B_2.standard),
                                                     0.0,accuracy: accuracy)
        
        
        let a2 = A_2.drop
        let b2 = B_2.drop
        
        XCTAssertEqual(a.x,a2.x,accuracy: accuracy)
        XCTAssertEqual(a.y,a2.y,accuracy: accuracy)
        XCTAssertEqual(b.x,b2.x,accuracy: accuracy)
        XCTAssertEqual(b.y,b2.y,accuracy: accuracy)

    }
     
     
}
