import { createFileRoute } from "@tanstack/react-router";
import Characterize from "../components/menu/characterize";

export const Route = createFileRoute("/characterize")({ component: Characterize });
