﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.MeshPartitioning
{
    public class GuidedPartioner2D<TElement> where TElement : class, ICell<INode>
    {
        private readonly IMesh2D<INode, TElement> mesh;
        private readonly Dictionary<int, IRegion2D> regions;

        public GuidedPartioner2D(IMesh2D<INode, TElement> mesh, Dictionary<int, IRegion2D> regions)
        {
            this.mesh = mesh;
            this.regions = regions;
        }

        public Dictionary<int, List<TElement>> CreateSubdomains()
        {
            var partition = new Dictionary<int, List<TElement>>();
            foreach (int id in regions.Keys) partition[id] = new List<TElement>();
            var partitionedElements = new Dictionary<TElement, int>(); //TODO: This should be replaced by accessing IElement.Subdomain
            var boundaryElements = new HashSet<TElement>();

            // Process elements that are unambiguously internal to one subdomain.
            for (int e = 0; e < mesh.Elements.Count; ++e)
            {
                TElement element = mesh.Elements[e];
                List<int> containingSubdomains = FindRegionsContainingElement(element);
                if (containingSubdomains.Count == 0)
                {
                    throw new Exception($"Element {e} does not belong to any of the provided regions.");
                }
                else if (containingSubdomains.Count == 1)
                {
                    int subdomain = containingSubdomains[0];
                    partition[subdomain].Add(element);
                    partitionedElements[element] = subdomain;
                }
                else
                {
                    boundaryElements.Add(element);
                }
            }

            // Process elements that are intersected by the region boundaries
            foreach (TElement element in boundaryElements)
            {
                int subdomain = FindSubdomainWithMostNeighbors(element, partitionedElements);
                partition[subdomain].Add(element);
                partitionedElements[element] = subdomain; // Take this element into account for the next boundary elements.
            }

            return partition;
        }

        private List<TElement> FindElementsWithCommonEdges(TElement element)
        {
            //TODO: At least 2 common nodes works for 1st order elements, but what about 2nd order ones?
            //TODO: Perhaps I should implement a method to get the Edge of an element and work with that one.

            // Find all elements with at least one common node. Also store the common nodes.
            var allNeighbors = new Dictionary<TElement, HashSet<INode>>();
            foreach (INode node in element.Nodes)
            {
                foreach (TElement neighbor in mesh.FindElementsWithNode(node))
                {
                    bool neighborIsKnown = allNeighbors.TryGetValue(neighbor, out HashSet<INode> commonNodes);
                    if (!neighborIsKnown) 
                    {
                        commonNodes = new HashSet<INode>();
                        allNeighbors.Add(neighbor, commonNodes);
                    }
                    commonNodes.Add(node);
                }
            }

            // Only keep the elements that have 2 or more common nodes.
            var elementsWithCommonEdge = new List<TElement>();
            foreach (var neighborNodesPair in allNeighbors)
            {
                if (neighborNodesPair.Value.Count > 1) elementsWithCommonEdge.Add(neighborNodesPair.Key);
            }

            return elementsWithCommonEdge;
        }

        private List<int> FindRegionsContainingElement(TElement element)
        {
            var containingRegionIds = new List<int>();
            foreach (var idRegionPair in regions)
            {
                int id = idRegionPair.Key;
                IRegion2D region = idRegionPair.Value;
                foreach (INode node in element.Nodes) //TODO: If a node was internal to a previous region, there is no need to process it again
                {
                    // We are interested in regions where nodes are internal to. If an element has 3 nodes internal to region A 
                    // and 1 node on the boundary between regions A,B then it belongs to A.
                    NodePosition relativePosition = region.FindRelativePosition(node);
                    if (relativePosition == NodePosition.Internal)
                    {
                        containingRegionIds.Add(id);
                        break;
                    }
                }
            }
            return containingRegionIds;
        }

        private int FindSubdomainWithMostNeighbors(TElement element, Dictionary<TElement, int> currentlyPartitionedElements)
        {
            // Count how many neighboring elements are contained in each subdomain 
            List<TElement> neighbors = FindElementsWithCommonEdges(element);
            var subdomainHistogram = new Dictionary<int, int>(); // Key = subdomain ID, Value = frequency
            foreach (int subdomain in regions.Keys) subdomainHistogram[subdomain] = 0;
            foreach (TElement neighbor in neighbors)
            {
                // The neighbor could be on the boundary. 
                bool isPartitioned = currentlyPartitionedElements.TryGetValue(element, out int subdomain);
                if (isPartitioned) subdomainHistogram[subdomain] += 1;
            }

            // Key = subdomain ID, Value = frequency
            IOrderedEnumerable<KeyValuePair<int, int>> sortedSubdomains =
                    subdomainHistogram.OrderByDescending(subdomainFrequencyPair => subdomainFrequencyPair.Value);
            Debug.Assert(sortedSubdomains.First().Value != 0, "None of this element's neighbors has been partitioned yet");
            return sortedSubdomains.First().Key;
        }
    }
}
