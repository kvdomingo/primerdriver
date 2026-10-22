import { createFileRoute } from "@tanstack/react-router";
import ProteinView from "../components/menu/protein-view";

export const Route = createFileRoute("/protein")({ component: ProteinView });
