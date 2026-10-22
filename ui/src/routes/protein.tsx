import { createFileRoute } from "@tanstack/react-router";
import { $api, queryClient } from "../api";
import ProteinView from "../components/menu/protein-view";

export const Route = createFileRoute("/protein")({
  loader: () =>
    queryClient.query({
      ...$api.queryOptions("get", "/api/expressionsys"),
      staleTime: "static",
    }),
  component: ProteinView,
});
