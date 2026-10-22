import { QueryClient } from "@tanstack/react-query";
import createFetchClient from "openapi-fetch";
import createClient from "openapi-react-query";
import type { components, paths } from "./openapi.gen";

export type ResultData =
  | components["schemas"]["DesignResult"]
  | components["schemas"]["MessageResponse"];

export const $api = createClient(createFetchClient<paths>());

export const queryClient = new QueryClient({
  defaultOptions: { queries: { staleTime: Number.POSITIVE_INFINITY } },
});
