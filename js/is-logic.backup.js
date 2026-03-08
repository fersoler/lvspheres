/**
 * Port of ISbuild.R and related logic to JavaScript
 */

function getSubPosition(sub) {
    return sub.reduce((acc, val) => acc + Math.pow(2, val - 1), 0) + 1;
}

function getSolPosition(sol) {
    let acc = 0;
    for (let i = 0; i < sol.length; i++) {
        if (sol[i] > 1e-10) {
            acc += Math.pow(2, i);
        }
    }
    return acc + 1;
}

function getAllSuperPos(list1, list2) {
    const l2 = list2.length;
    const v2 = getSubPosition(list1) - 1;
    let results = [];
    for (let x = 1; x < Math.pow(2, l2); x++) {
        let prov = [];
        for (let i = 0; i < l2; i++) {
            if ((x >> i) & 1) prov.push(list2[i]);
        }
        results.push(getSubPosition(prov) + v2);
    }
    return results;
}

function getCombinations(array, size) {
    function p(t, a, s) {
        if (t.length === s) {
            results.push(t);
            return;
        }
        for (let i = 0; i < a.length; i++) {
            p(t.concat(a[i]), a.slice(i + 1), s);
        }
    }
    const results = [];
    p([], array, size);
    return results;
}

function ISbuild(alphas0, gammas) {
    const size = alphas0.length;
    let allStatPoints = [];
    let foundSols = 0;

    // Trivial solution
    allStatPoints.push(new Array(size).fill(0));
    foundSols++;

    let listSubsGASS = new Array(Math.pow(2, size)).fill(0);
    listSubsGASS[0] = 1;

    let listGASSindex = new Array(Math.pow(2, size)).fill(0);
    listGASSindex[0] = 1;

    let subsetGASS = [];

    for (let s = size; s >= 1; s--) {
        const combinations = getCombinations(Array.from({ length: size }, (_, i) => i + 1), s);
        for (let i = combinations.length - 1; i >= 0; i--) {
            const ss = combinations[i];
            const subsCode = getSubPosition(ss);

            if (listGASSindex[subsCode - 1] === 0) {
                // Solve LCP: q = -alphas, M = -gammas (Wait, R uses matrix(-gammas) and q = -alphas)
                const subAlphas = ss.map(sIdx => alphas0[sIdx - 1]);
                const subGammas = ss.map(rIdx => ss.map(cIdx => gammas[rIdx - 1][cIdx - 1]));

                // Lemke expects w - Mz = q. 
                // For biodiversity: gammas*z + alphas = w  => w - gammas*z = alphas
                // So M = gammas, q = alphas. (Check lemkelcp.R: matrix(-gammas), t(-alphas))
                // R lemkelcp(M, q) solves w - Mz = q.
                // In R: lemkelcp(matrix(-gammas), t(-alphas)) => w - (-gammas)z = -alphas => w + gammas*z = -alphas => gammas*z + alphas = -w
                // Actually LV is: diag(z)(alphas + gammas*z) = 0.
                // Equilibrium: alphas + gammas*z <= 0, and if z_i > 0, (alphas + gammas*z)_i = 0.
                // This is LCP(M, q) with M = -gammas, q = -alphas.
                // w = (-gammas)z + (-alphas) >= 0, z >= 0, w^T z = 0.
                const M_lcp = subGammas.map(row => row.map(v => -v));
                const q_lcp = subAlphas.map(v => -v);

                const sol = lemkelcp(M_lcp, q_lcp);
                if (sol.z) {
                    const gass = sol.z;
                    const v = new Array(size).fill(0);
                    ss.forEach((sIdx, k) => v[sIdx - 1] = gass[k]);
                    const gassCode = getSolPosition(v);

                    let lookAt = listSubsGASS[gassCode - 1];
                    if (lookAt === 0) {
                        foundSols++;
                        lookAt = foundSols;
                        allStatPoints.push(v);
                        listSubsGASS[gassCode - 1] = lookAt;
                    }
                    listGASSindex[subsCode - 1] = lookAt;
                    subsetGASS.push({ subset: ss.join(", "), ind: lookAt });

                    // The "look for sub-communities" logic is an optimization for speed in R,
                    // we can skip it for a 3-species system in JS as it's almost instant anyway.
                }
            }
        }
    }

    // Step 2: Connectivity
    let connectivity = Array.from({ length: foundSols }, () => new Array(foundSols).fill(0));
    for (let p = 0; p < foundSols; p++) {
        const point = allStatPoints[p];
        let I = [];
        let NI = [];
        for (let i = 0; i < size; i++) {
            if (point[i] > 1e-10) I.push(i + 1);
            else NI.push(i + 1);
        }

        if (NI.length > 0) {
            let R = [];
            if (I.length === 0) {
                R = NI.map(idx => alphas0[idx - 1]);
            } else {
                // R = alphas[NI] + gammas[NI,I] %* % point[I]
                R = NI.map(ri => {
                    let sum = alphas0[ri - 1];
                    for (let j = 0; j < I.length; j++) {
                        const ci = I[j];
                        sum += gammas[ri - 1][ci - 1] * point[ci - 1];
                    }
                    return sum;
                });
            }

            let J = [];
            for (let k = 0; k < R.length; k++) {
                if (R[k] > 1e-10) J.push(NI[k]);
            }

            if (J.length > 0) {
                const superPos = getAllSuperPos(I, J);
                superPos.forEach(spCode => {
                    const targetIdx = listGASSindex[spCode - 1];
                    if (targetIdx > 0) connectivity[p][targetIdx - 1] = 1;
                });
            }
        }
    }

    const gassInd = listGASSindex[getSubPosition(Array.from({ length: size }, (_, i) => i + 1)) - 1];

    return {
        points: allStatPoints,
        connectivity: connectivity,
        gassInd: gassInd
    };
}

// Layout logic
const vertexCoords3D = {
    "0,0,0": [0, 0],
    "1,0,0": [2, 0],
    "0,1,0": [0, 2],
    "1,1,0": [2, 2],
    "0,0,1": [-1, -1.5],
    "1,0,1": [1, -1.5],
    "0,1,1": [-1, 0.5],
    "1,1,1": [1, 0.5]
};

function getCoordIS3D(points) {
    return points.map(p => {
        const key = p.map(v => (v > 1e-10 ? 1 : 0)).join(",");
        return vertexCoords3D[key] || [0, 0];
    });
}
