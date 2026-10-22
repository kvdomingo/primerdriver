import { createFileRoute } from "@tanstack/react-router";
import DnaView from "../components/menu/dna-view";

export const Route = createFileRoute("/dna")({ component: DnaView });
